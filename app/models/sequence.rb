class Sequence 
  include DataMapper::Resource
  
  property :seq_id, Serial
  property :seq_name, String, :required => true
  property :sequence, Text, :required => false
  property :seq_type, String, :required => true
  property :seq_accession, String, :required => true
  property :abrev_name, String, :required => false
  property :disorder_percent, Integer, :required => false
  
  #has n, :a_asequence
  #has n, :disorder
end

