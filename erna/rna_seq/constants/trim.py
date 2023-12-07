TRIM = {
  # minimum length after trimming
  'min_len': 12,
  # adapter sequences
  'adapters': {
    'NEB': {
      'adapter_3end': 'AGATCGGAAGAGCACACGTCT',
      'adapter_5end': 'GTTCAGAGTTCTACAGTCCGACGATC',    
    },
    'Illumina': {
      'adapter_5end': 'GTTCAGAGTTCTACAGTCCGACGATC',
      'adapter_3end': 'TGGAATTCTCGGGTGCCAAGG', 
    },
    # TODO: confirm the adapter name
    'QG': {
      'adapter_3end': 'AACTGTAGGCACCATCAAT',
    },
  },
  # trim sequence by fixed length
  'trim': {
    'keep_5end': 12,
    'keep_3end': 12,
  },
}