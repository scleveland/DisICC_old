class Disorder
  include DataMapper::Resource
  
  property :disorder_id, Serial
  property :seq_id, Integer, :required => true
  property :disorder_type, String, :required => true
  
  
end
