

class Sequence 
  include DataMapper::Resource
  
  property :seq_id, Serial
  property :seq_name, String, :required => true
  property :sequence, Text, :required => false
  property :seq_type, String, :required => true
  property :seq_accession, String, :required => true
  property :abrev_name, String, :required => false
  property :disorder_percent, Integer, :required => false
  property :alternate_name, String, :required => false
  
  
  #has n, :a_asequences
  #has n, :disorder
  
  def generate_fasta_file
    filepath = "temp_data/"+self.abrev_name+"_"+self.seq_type+".fasta"
    f = File.new(filepath, "w+")
    f.write(">"+self.abrev_name+"|"+self.seq_name+"|"+self.seq_type+"|"+self.seq_accession+"\n")
    f.write(self.sequence)
    f.close
    return filepath
  end
  
  def generate_iupred_disorder
    res = `./lib/disorder_apps/iupred/iupred_mac #{self.generate_fasta_file} long`
    filepath = "temp_data/"+self.abrev_name+"_"+self.seq_type+"_iupred.fasta"
    f = File.new(filepath, "w+")
    f.write(res)
    f.close
    self.store_iupred(filepath)
  end
  
  def store_iupred(filepath)
    #create a new disorder object
    dis = Disorder.create(:seq_id => self.seq_id, :disorder_type=>"IUPred", :version=>1)
    file = File.new(filepath, 'r')
    counter = 1
    aa_count = 0
    while (line = file.gets)
      #puts "#{counter}: #{line}"
      counter = counter + 1
      if counter > 10
        line_array = line.split(' ')
        if aa = AAsequence.first(:seq_id => self.seq_id, :original_position=>aa_count, :amino_acid=>line_array[1])
        #puts "Amino Acid -"+line_array[1]+ " : " + aa.amino_acid + " | " + aa_count.to_s
        #if aa.amino_acid == line_array[1]  
          dv = DisorderValue.create(:disorder_id => dis.disorder_id, :aasequence_id => aa.AAsequence_id, :dvalue=>line_array[2].to_f) 
        end
        aa_count +=1
      end
    end
  end
  
  
  # jalview_string = ""
  # low disorder  ffff00        
  # avg disorder  ffcc00        
  # medium disorder ff9900        
  # highly disordered ff6600        
  # extremely disordered  ff0000        
  # no disorder 0
  
  def generate_all_sequence_javliew_annotation_iupred(ptype)
    require 'csv'
    CSV.open("temp_data/jalview_#{ptype}_#{Time.now}.gff", "wb", {:col_sep => "\t"}) do |csv|
    csv <<["possible disorder","ffffcc"]
    csv <<["low disorder","ffff00"]
    csv <<["avg disorder","ffcc00"]
    csv <<["medium disorder","ff9900"]
    csv <<["highly disordered","ff6600"]
    csv <<["extremely disordered","ff0000"]
    csv <<["no disorder","0"]
    Sequence.all(:seq_type => ptype).each do |seq|
      if dis = Disorder.first(:seq_id=>seq.seq_id, :disorder_type=>"IUPred")
      counter = 1
      DisorderValue.all(:disorder_id=>dis.disorder_id,:dvalue.gte => 0.4, :order=>[:disorder_value_id]).each do |dv|
        if dv.dvalue > 0.4 && dv.dvalue < 0.5
          feature_type = "possible disorder"
        elsif dv.dvalue > 0.5 && dv.dvalue < 0.6
          feature_type = "low disorder"
        elsif dv.dvalue >= 0.6 && dv.dvalue < 0.7
          feature_type = "avg disorder"
        elsif dv.dvalue >= 0.7 && dv.dvalue < 0.8
          feature_type = "medium disorder"
        elsif dv.dvalue >= 0.8 && dv.dvalue < 0.9
          feature_type = "highly disordered"
        elsif dv.dvalue > 0.9
          feature_type = "extremely disordered"
        else
          feature_type = "no disorder"
        end
        unless seq.alternate_name.nil?
          
          csv << ["None", "#{seq.alternate_name.split('/')[0]}", -1, "#{dv.aasequence.original_position}", "#{dv.aasequence.original_position}", "#{feature_type}"]
        end
        counter+=1
      end
      end
     end
    end
  end
  
  def generate_jalview_annotation_iupred
    jalview_string= ""
    dis = Disorder.first(:seq_id=>self.seq_id)
    counter = 1
    DisorderValue.all(:disorder_id=>dis.disorder_id,:dvalue.gte => 0.5, :order=>[:disorder_value_id]).each do |dv|
      if dv.dvalue > 0.5 && dv.dvalue < 0.6
        feature_type = "low disorder"
      elsif dv.dvalue >= 0.6 && dv.dvalue < 0.7
        feature_type = "avg disorder"
      elsif dv.dvalue >= 0.7 && dv.dvalue < 0.8
        feature_type = "medium disorder"
      elsif dv.dvalue >= 0.8 && dv.dvalue < 0.9
        feature_type = "highly disordered"
      elsif dv.dvalue > 0.9
        feature_type = "extremely disordered"
      else
        feature_type = "no disorder"
      end
      jalview_string = jalview_string + "None #{self.abrev_name} -1 #{counter} #{counter} #{feature_type}"+"\n"
   end
   return jalview_string
  end
end

