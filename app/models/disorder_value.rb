class DisorderValue
  include DataMapper::Resource
  
  property :disorder_value_id, Serial
  property :disorder_id, Integer, :required => true
  property :aasequence_id, Integer, :required => true
  property :dvalue, Float, :required => true

  

end
