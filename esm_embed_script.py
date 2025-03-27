#!/usr/bin/env python3

import os
import torch  # type: ignore
import esm  # type: ignore
import argparse
import numpy as np  # type: ignore
import pandas as pd  # type: ignore
from tqdm import tqdm  # type: ignore
from typing import List, Dict, Tuple, Optional


def calculate_entropy(counts: pd.Series) -> float:
    """
    Calculate Shannon entropy (in bits) from a series of counts.

    Args:
        counts: Pandas Series of sequence counts

    Returns:
        Shannon entropy using base 2 logarithm.
    """
    total = counts.sum()
    if total == 0:
        return 0.0
    proportions = counts / total
    non_zero = proportions > 0
    return -np.sum(proportions[non_zero] * np.log2(proportions[non_zero]))


def load_esm2_model(model_name: str = "esm2_t33_650M_UR50D") -> Tuple:
    """
    Load an ESM-2 model and tokenizer
    
    Args:
        model_name: Name of the ESM-2 model to use
        
    Returns:
        Tuple of (model, alphabet, batch_converter)
    """
    print(f"Loading ESM-2 model: {model_name}")
    
    # Load the model
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model.eval()  # Set to evaluation mode
    
    # Check if GPU is available and move model to GPU
    if torch.cuda.is_available():
        model = model.cuda()
        print("Using GPU for inference")
    else:
        print("Using CPU for inference")
    
    # Create batch converter
    batch_converter = alphabet.get_batch_converter()
    
    return model, alphabet, batch_converter


def load_metadata(metadata_file: str) -> pd.DataFrame:
    """
    Load patient metadata
    
    Args:
        metadata_file: Path to metadata TSV file
        
    Returns:
        DataFrame with patient metadata
    """
    # Load metadata
    metadata = pd.read_csv(metadata_file, sep='\t', header=0)
    
    # Rename columns for clarity
    metadata.columns = ['patient_id', 't1d_status', "batch_name"]
    metadata.patient_id = metadata.patient_id.str.removesuffix(
        '.tsv' if any( metadata.patient_id.map(lambda x: '.tsv' in x) ) else '.csv'
        )
    
    n_t1d = sum([1 for i in metadata['t1d_status'] if i == 'Yes'])
    print(f"Loaded metadata for {len(metadata)} patients")
    print(f"T1D patients: {n_t1d}")
    print(f"Healthy controls: {len(metadata) - n_t1d}")
    
    return metadata


def load_patient_repertoire(repertoire_file: str, count_threshold: Optional[int] = None) -> pd.DataFrame:
    """
    Load TCR repertoire for a single patient
    
    Args:
        repertoire_file: Path to patient repertoire TSV file
        count_threshold: Minimum number of times a clone must appear in the repertoire to be included in the dataframe
        
    Returns:
        DataFrame with TCR sequences and counts
    """
    # Load repertoire data
    repertoire = pd.read_csv(repertoire_file, sep='\t', header=0)
    
    # Rename columns for clarity
    repertoire.columns = ['cdr3_aa', 'n_seq']

    # Merge and sum the same cdr3
    repertoire = repertoire.groupby('cdr3_aa', as_index=False).\
        sum().\
        sort_values(by='n_seq', ascending=False).\
        reset_index(drop=True)
    
    # Extract patient ID from filename
    patient_id = os.path.basename(repertoire_file).split('.')[0]
    
    # Add patient ID column
    repertoire['patient_id'] = patient_id

    # Filter out clones with small number of counts
    if count_threshold:
        repertoire = repertoire[repertoire['n_seq'] >= count_threshold]
    
    return repertoire


def generate_embeddings_batch(model, batch_converter, sequences: List[str], 
                             ids: List[str], batch_size: int = 32, 
                             repr_layer: int = 33) -> Dict:
    """
    Generate embeddings for a batch of sequences
    
    Args:
        model: ESM-2 model
        batch_converter: Batch converter for the model
        sequences: List of amino acid sequences
        ids: List of sequence identifiers
        batch_size: Batch size for processing
        repr_layer: Which layer to extract representations from
        
    Returns:
        Dictionary mapping sequence IDs to embeddings
    """
    results = {}
    
    # Process in batches
    for i in range(0, len(sequences), batch_size):
        batch_ids = ids[i:i+batch_size]
        batch_seqs = sequences[i:i+batch_size]
        
        # Prepare batch data
        batch_labels, batch_strs, batch_tokens = batch_converter(list(zip(batch_ids, batch_seqs)))
        
        # Move to GPU if available
        if torch.cuda.is_available():
            batch_tokens = batch_tokens.cuda()
        
        with torch.no_grad():
            # Extract embeddings from the model
            token_representations = model(batch_tokens, repr_layers=[repr_layer])["representations"][repr_layer]
        
        # Store results for each sequence
        for j, seq_id in enumerate(batch_ids):
            # Get sequence length (exclude start/end tokens)
            seq_len = len(batch_seqs[j])
            
            # Get representations (exclude start/end tokens, padding)
            seq_repr = token_representations[j, 1:seq_len+1].cpu().numpy()
            
            # Store in results dictionary
            results[seq_id] = seq_repr
    
    return results


def compute_sequence_embedding(token_embeddings: Dict, pool_method: str = "mean") -> Dict:
    """
    Compute sequence-level embeddings from token-level embeddings
    
    Args:
        token_embeddings: Dictionary with token-level embeddings
        pool_method: Method for pooling ('mean' or 'cls')
        
    Returns:
        Dictionary with sequence-level embeddings
    """
    sequence_embeddings = {}
    
    for seq_id, token_emb in token_embeddings.items():
        if pool_method == "mean":
            # Mean pooling across all tokens
            sequence_embeddings[seq_id] = np.mean(token_emb, axis=0)
        elif pool_method == "cls":
            # Use first token as representation
            sequence_embeddings[seq_id] = token_emb[0]
        else:
            raise ValueError(f"Unknown pooling method: {pool_method}")
    
    return sequence_embeddings


def process_patient_repertoire(patient_id: str, repertoire_df: pd.DataFrame, 
                              model, alphabet, batch_converter,
                              batch_size: int = 32, 
                              pool_method: str = "mean",
                              repr_layer: int = 33) -> np.ndarray:
    """
    Process TCR repertoire for a single patient and generate weighted embedding
    
    Args:
        patient_id: Patient identifier
        repertoire_df: DataFrame with patient's TCR repertoire
        model: ESM-2 model
        alphabet: ESM-2 alphabet
        batch_converter: Batch converter for the model
        batch_size: Batch size for processing
        pool_method: Method for pooling token embeddings
        repr_layer: Which layer to extract representations from
    
    Returns:
        Weighted average embedding for the patient
    """
    # Extract sequences and their counts
    sequences = repertoire_df['cdr3_aa'].tolist()
    counts = repertoire_df['n_seq'].tolist()
    
    # Create unique sequence IDs
    seq_ids = [f"{patient_id}_{i}" for i in range(len(sequences))]
    
    # Filter out invalid sequences
    valid_indices = []
    valid_sequences = []
    valid_ids = []
    valid_counts = []
    
    for i, (seq_id, seq, count) in enumerate(zip(seq_ids, sequences, counts)):
        if pd.notna(seq) and all(aa in alphabet.all_toks for aa in seq):
            valid_indices.append(i)
            valid_sequences.append(seq)
            valid_ids.append(seq_id)
            valid_counts.append(count)
    
    if not valid_sequences:
        print(f"Warning: No valid sequences found for patient {patient_id}")
        # Return zero vector of appropriate dimension
        if model.embed_dim:
            return np.zeros(model.embed_dim)
        else:
            return np.zeros(1280)  # Default for esm2_t33_650M_UR50D
    
    print(f"Processing {len(valid_sequences)} valid sequences for patient {patient_id}")
    
    # Generate token-level embeddings
    token_embeddings = generate_embeddings_batch(
        model, batch_converter, valid_sequences, valid_ids, 
        batch_size=batch_size, repr_layer=repr_layer
    )
    
    # Compute sequence-level embeddings
    sequence_embeddings = compute_sequence_embedding(token_embeddings, pool_method)
    
    # Convert to numpy arrays for weighted average
    embedding_arrays = []
    weights = []
    
    for i, seq_id in enumerate(valid_ids):
        embedding_arrays.append(sequence_embeddings[seq_id])
        weights.append(valid_counts[i])
    
    embedding_matrix = np.stack(embedding_arrays)
    weight_array = np.array(weights)
    
    # Normalize weights
    normalized_weights = weight_array / weight_array.sum()
    
    # Compute weighted average
    weighted_embedding = np.sum(embedding_matrix * normalized_weights[:, np.newaxis], axis=0)
    
    return weighted_embedding


def process_all_patients(data_dir: str, output_file: str, 
                         model_name: str = "esm2_t33_650M_UR50D",
                         batch_size: int = 32, pool_method: str = "mean",
                         repr_layer: int = 33):
    """
    Process TCR repertoires for all patients and generate embeddings
    
    Args:
        data_dir: Directory containing metadata and repertoire files
        output_file: File name to save patients' embeddings into without file extension spec
        model_name: ESM-2 model to use
        batch_size: Batch size for processing
        pool_method: Method for pooling token embeddings
        repr_layer: Which layer to extract representations from
    """
    # Load metadata
    metadata_file = os.path.join(data_dir, "metadata.tsv")
    metadata = load_metadata(metadata_file)
    
    # Load ESM-2 model
    model, alphabet, batch_converter = load_esm2_model(model_name)
    
    # Process each patient
    patient_embeddings = []
    repertoire_dir = os.path.join(data_dir, "repertoires")
    
    for i, patient_id in enumerate(tqdm(metadata['patient_id'], desc="Processing patients")):
        # Load patient repertoire
        repertoire_file = os.path.join(repertoire_dir, f"{patient_id}.tsv")
        
        if not os.path.exists(repertoire_file):
            print(f"Warning: Repertoire file not found for patient {patient_id}")
            continue
        
        repertoire_df = load_patient_repertoire(repertoire_file)

        # Calculate entropy
        entropy = calculate_entropy(repertoire_df['n_seq'])
        
        # Generate embedding
        patient_embedding = process_patient_repertoire(
            patient_id, repertoire_df, model, alphabet, batch_converter,
            batch_size=batch_size, pool_method=pool_method, repr_layer=repr_layer
        )
        
        # Add to results
        t1d_status = metadata.loc[metadata['patient_id'] == patient_id, 't1d_status'].iloc[0]
        batch_name = metadata.loc[metadata['patient_id'] == patient_id, 'batch_name'].iloc[0]
        
        patient_embeddings.append({
            'patient_id': patient_id,
            't1d_status': t1d_status,
            'batch_name': batch_name,
            'embedding': patient_embedding,
            'num_sequences': len(repertoire_df),
            'entropy': entropy,
        })
        
        # Save intermediate results periodically
        if (i + 1) % 10 == 0 or (i + 1) == len(metadata):
            temp_df = pd.DataFrame(patient_embeddings)
            temp_file = os.path.join(data_dir, f"{output_file}_temp.pkl")
            temp_df.to_pickle(temp_file)
            print(f"Saved intermediate results to {temp_file}")
    
    # Create final DataFrame
    result_df = pd.DataFrame(patient_embeddings)
    
    # Save results
    result_file = os.path.join(data_dir, f"{output_file}.pkl")
    result_df.to_pickle(result_file)
    print(f"Saved patient embeddings to {output_file}")
    
    # Add embedding dimensions
    temp_dict = {}
    embedding_dim = model.embed_dim
    for i in range(embedding_dim):
        temp_dict[f'dim_{i}'] = result_df['embedding'].apply(lambda x: x[i])

    # Make resulting DataFrame
    result_csv = pd.DataFrame({
        'patient_id': result_df['patient_id'],
        't1d_status': result_df['t1d_status'],
        'batch_name': result_df['batch_name'],
        'num_sequences': result_df['num_sequences'],
        'entropy': result_df['entropy'],
        **temp_dict
    })
    
    csv_output = os.path.join(data_dir, f"{output_file}.csv")
    result_csv.to_csv(csv_output, index=False)
    print(f"Saved patient embeddings CSV to {csv_output}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate patient-level TCR embeddings from repertoire files"
    )
    parser.add_argument(
        "--data-dir",
        required=True,
        help="Directory containing metadata.tsv and repertoire folder",
    )
    parser.add_argument(
        "--output", required=True, help="Output file for patient embeddings"
    )
    parser.add_argument(
        "--model", default="esm2_t33_650M_UR50D", help="ESM-2 model to use"
    )
    parser.add_argument(
        "--batch-size", type=int, default=32, help="Batch size for processing"
    )
    parser.add_argument(
        "--pool-method",
        default="mean",
        choices=["mean", "cls"],
        help="Method for pooling token embeddings",
    )
    parser.add_argument(
        "--repr-layer",
        type=int,
        default=33,
        help="Model layer to extract representations from",
    )

    args = parser.parse_args()

    # Process all patients
    process_all_patients(
        data_dir=args.data_dir,
        output_file=args.output,
        model_name=args.model,
        batch_size=args.batch_size,
        pool_method=args.pool_method,
        repr_layer=args.repr_layer,
    )


if __name__ == "__main__":
    main()  # pragma: no cover
