class Alignment
  include DataMapper::Resource
  
  property :align_id, Serial
  property :seq_id, Integer, :required => true
  property :alignment_name, String, :required => true
  property :align_order, Integer, :required => true
  property :alignment_sequence, Text, :required => true, :default => ""
  # has n, :disorder_values
  # belongs_to :sequence
  
end
