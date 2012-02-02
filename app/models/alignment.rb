class Alignment
  include DataMapper::Resource
  
  property :align_id, Serial
  property :seq_id, Integer, :required => true
  property :alignment_name, String, :required => true
  property :align_order, Integer, :required => true
  property :alignment_sequence, Text, :required => true, :default => ""
  property :fasta_title, Text, :required => false, :default => ""
  # has n, :disorder_values
  # belongs_to :sequence
  
  validates_uniqueness_of :alignment_name
  
  #runs the score app to get percent identities for each sequence
  def run_score
    filename = self.generate_fasta_alignment_file
    string = "./lib/score_mac #{filename} temp_data/#{self.alignment_name}_res.txt temp_data/#{self.alignment_name}_dif.txt temp_data/#{self.alignment_name}_alignments.txt"
    puts string
    if system(string)
      
    end
  end
  
  def sequences
    alignments = Alignment.all(:alignment_name => self.alignment_name, :order=>[:align_order])
    seq_ids = alignments.map{|a| a.seq_id}
    Sequence.all(:seq_id=>seq_ids)
  end
  
  
  #runs the AlignAssess app to get percent identities for each sequence in the alignment
  def run_align_assess
    filename = self.generate_fasta_alignment_file_for_all
    string = "./lib/AlignAssess_wShorterID #{filename} P"
    seq_array = Array.new
    if system(string)
      seq_id_array = self.sequences.map{|s| s.seq_id}
      new_filename = filename + "_assess"
      f = File.new(new_filename, "r")
      flag = false
      read_row= 999999999
      cur_row = 0
      while (line = f.gets)
        if cur_row > read_row  && flag
          if line == "\n"
            flag =false
          else
            seq_array << line.split("\t")
          end
        elsif line == "Pair-wise %ID over shorter sequence:\n"
          flag=true
          read_row = cur_row + 2
        end
        cur_row +=1
      end
      range = seq_array.length - 1
      #seq_array.each do |row|
      for row_num in 0..range
        for i in 1..range#(row_num)          
          PercentIdentity.first_or_create(:seq1_id=>seq_id_array[row_num],
                                                    :seq2_id=>seq_id_array[i],
                                                    :alignment_name => self.alignment_name,
                                                    :percent_id=>seq_array[row_num][i])
          # print "[#{row_num}:#{i-1}=>#{row[i]}],"
        end
        #print "\n"
      end
    end
  end
  
  def run_caps
    self.run_align_assess
    Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    alignments.each do |alignment|
      alignment.generate_pid_fasta_file("temp_data/#{self.alignment_name}")
    end
    system "./lib/comp_apps/caps/caps -F temp_data/#{self.alignment_name} --intra"
  end
  
  def run_xdet
    self.run_align_assess
    Dir.mkdir("temp_data/#{self.alignment_name}") unless File.directory?("temp_data/#{self.alignment_name}")
    alignments = Alignment.all(:alignment_name => self.alignment_name)
    alignments.each do |alignment|
      filename= alignment.generate_pid_fasta_file("temp_data/#{self.alignment_name}")
      system "./lib/comp_apps/XDet/xdet_linux32 #{filename} ~/Rails/DisICC/lib/comp_apps/XDet/Maxhom_McLachlan.metric >> #{filename}_xdet"
    end
  end
  
  def generate_pid_fasta_files
    self.run_align_assess
    self.sequences.each do |seq|
      fasta_string=""
      pids = PercentIdentity.all(:seq1_id => seq.seq_id, :percent_id.gte => 20, :order =>[:percent_id.desc])
      fasta_string= Alignment.first(:alignment_name => self.alignment_name, :seq_id=>seq.seq_id).fasta_alignment_string
      puts seq.abrev_name+":"+pids.count.to_s
      pids.each do |pid|
        if pid.seq2_id != seq.seq_id
          print Sequence.get(pid.seq2_id).abrev_name + ":" + pid.percent_id.to_s + ","
          fasta_string = fasta_string + Alignment.first(:alignment_name=>pid.alignment_name, :seq_id=>pid.seq2_id).fasta_alignment_string("pid:#{pid.percent_id}")
        end
      end
      filepath = "temp_data/"+self.alignment_name+"_"+seq.abrev_name+"_pid.fasta"
      f = File.new(filepath, "w+")
      f.write(fasta_string)
      f.close
      print "\n"
    end
  end
  
  #assumes that run_align_asses has already occured
  def generate_pid_fasta_file(dir="temp_data")
    fasta_string=""
    seq = Sequence.get(self.seq_id)
    pids = PercentIdentity.all(:seq1_id => self.seq_id, :percent_id.gte => 20, :order =>[:percent_id.desc])
    fasta_string= Alignment.first(:alignment_name => self.alignment_name, :seq_id=>self.seq_id).fasta_alignment_string
    puts seq.abrev_name+":"+pids.count.to_s
    pids.each do |pid|
      if pid.seq2_id != seq.seq_id
        print Sequence.get(pid.seq2_id).abrev_name + ":" + pid.percent_id.to_s + ","
        fasta_string = fasta_string + Alignment.first(:alignment_name=>pid.alignment_name, :seq_id=>pid.seq2_id).fasta_alignment_string("pid:#{pid.percent_id}")
      end
    end
    filepath = "#{dir}/"+self.alignment_name+"_"+seq.abrev_name+"_pid.fasta"
    f = File.new(filepath, "w+")
    f.write(fasta_string)
    f.close
    filepath
  end
  
  def update_alignment_sequence
    cur_position = 0   
    fasta_string = ""
    AlignmentPosition.all(:alignment_id => self.align_id, :order => [:alignmnet_position_id.asc]).each do |position|
      if position.position == cur_position
         fasta_string = fasta_string + AAsequence.first(:AAsequence_id => position.aasequence_id).amino_acid
      else
         while position.position > cur_position
                      fasta_string = fasta_string +"-"
                      cur_position += 1
         end
        fasta_string = fasta_string + AAsequence.first(:AAsequence_id => position.aasequence_id).amino_acid
      end
      cur_position += 1         
    end
    self.alignment_seqence = fasta_string
    self.save
  end
  
  def fasta_alignment_string(extra_string="")
    if self.alignment_sequence.empty?
      self.update_alignment_sequence
    end
    seq = Sequence.get(self.seq_id)
    fasta_string = ">"+seq.abrev_name+"|"+seq.seq_name+"|"+seq.seq_type+"|"+seq.seq_accession+"|"+extra_string+"\n" + self.alignment_sequence+"\n"
  end
  
  def generate_fasta_alignment_file_for_all
    alignments = Alignment.all(:alignment_name=> self.alignment_name)
    filepath = "temp_data/"+self.alignment_name+"_alignment"+Time.now.to_i.to_s+".fasta"
    f = File.new(filepath, "w+")
    alignments.each do |alignment|
      f.write(alignment.fasta_alignment_string)
    end
    f.close
    filepath
  end
  
 
  #This generates a fasta file for each sequence in the alignment
  #only sequences greater than 20% identitiy will be included
  def generate_fasta_files_for_comp_mut
  
  end
end
