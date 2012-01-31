class PercentIdentity
  include DataMapper::Resource
  
  property :seq1_id,  Integer, :key => true
  property :seq2_id, Integer, :key => true
  property :alignment_name, String, :key => true
  property :percent_id, Float, :required=> true
  
end