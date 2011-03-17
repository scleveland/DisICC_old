class AAsequence 
  include DataMapper::Resource
  
  property :AAsequence_id, Serial
  property :seq_id, Integer, :required => true
  property :amino_acid, String, :length=> 1,  :required => true
  property :original_position, Integer, :required => false
  property :disorder_consensus, Float, :required => false
  property :contact_consensus, Float, :required => false
  property :contact_positive_consensus, Integer, :required => false
  
  
end
