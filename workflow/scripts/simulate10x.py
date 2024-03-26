import numpy as np

#Parameter definition
large_droplet=40000000
class Molecule(object):
    def __init__(self,length,start,end,index):
        self.seqidx=index
        self.length=length
        self.index_droplet=0
        self.barcode='Null'
        self.shortreads=0
        self.start=start
        self.end=end

class Qual_Substitute(object):
    def __init__(self,phred,change):
        self.phred=phred
        self.substitute=change

class Short_reads_PE(object):
     def __init__(self,seq1,qual1,start1,end1,seq2,qual2,start2,end2):
         self.seq1=seq1
         self.qual1=qual1
         self.start1=start1
         self.end1=end1
         self.seq2=seq2
         self.qual2=qual2
         self.start2=start2
         self.end2=end2

class parameter(object):
    def __init__(self):
        self.fastaHap1='N'
        self.fastaHap2='N'
        self.barcodePool='N'
        self.barcodeQualityFile='N'
        self.seqQualityFile='N'
        self.coverageLongFrag=0
        self.avgLenLongFrag=0
        self.molPerDroplet=0
        self.lenShortRead=0
        self.coverageShortRead=0
        self.avgInsertShortRead=0
        self.stdInsertShortRead=0
        self.errorRate=0
        self.hap=1
        self.threads=1

#read template sequence
def input_seq(in_path):
    reflist = []
    reftitle = []
    sequence = ""
    compressed = True if in_path.lower().endswith("gz") else False
    if compressed:
        f = gzip.open(in_path, "r")
    else:
        f = open(in_path, "r")
    index = 0
    for line in f:
        line = line.decode().strip("\n") if compressed else line.strip("\n")
        if line[0] == ">":
            title = line.split(" ")[0]
            reftitle.append(title[1:])
            if index == 0:
                index += 1
            else:
                reflist.append(sequence)
                sequence = ""
        else:
            sequence += line
    reflist.append(sequence)
    f.close()
    return reflist, reftitle

def randomlong(Par,reflist,reftitle):
    """randomly generate long molecules from supplied genome""" 
    global N_frag
    N_frag=0
    frag_list=[]
    MolSet=[]
    For_hist=[]
    index=0
    T_frag=0
    T_frag_eff=0
    ref_total_len=0
    ref_effective_len=0
    for seq in reflist:
        #calculate required number of molecules
        lensingle=len(seq)
        N_frag=int(lensingle*Par.coverageLongFrag/(Par.avgLenLongFrag*1000))
        #randomly simulate molecules
        for i in range(N_frag):
            start=int(np.random.uniform(low=0,high=lensingle))
            length=int(np.random.exponential(scale=Par.avgLenLongFrag*1000))
            end=start+length-1
            if length==0:
               continue
            if end>lensingle:
               Molseq=seq[start:lensingle]
               lengthnew=lensingle-start
               NewMol=Molecule(lengthnew,start,lensingle,index)
               MolSet.append(NewMol)
            else:
               Molseq=seq[start:end]
               NewMol=Molecule(length-1,start,end,index)
               MolSet.append(NewMol)
        index+=1
    N_frag=len(MolSet)
    return MolSet

def deternumdroplet(N_frag,N_FP):
    """assign long fragments to droplet"""
    frag_drop = np.random.poisson(N_FP, large_droplet)
    assign_drop=[]
    totalfrag=0
    for i in range(large_droplet):
        totalfrag += frag_drop[i]
        if totalfrag<=N_frag:
            assign_drop.append(frag_drop[i])
        else:
            last=N_frag-(totalfrag-frag_drop[i])
            assign_drop.append(last)
            break
    return assign_drop

def selectbarcode(pool,assign_drop,MolSet,droplet_container):
    """assign barcode to each fragment"""
    #permute index of long fragment for random sampling
    permutnum=np.random.permutation(N_frag)
    #include barcode in the list
    barcode=[]
    N_droplet=len(assign_drop)
    t=1
    f=open(pool,"r")
    for line in f:
        barcode.append(line.strip('\n'))
        t += 1
        if t>N_droplet:
            break
    f.close()
    start=0
    for i in range(N_droplet):
        num_molecule_per_partition=assign_drop[i]
        index_molecule=permutnum[start:start+num_molecule_per_partition]
        totalseqlen=0
        temp=[]
        start += num_molecule_per_partition
        for j in range(num_molecule_per_partition):
            index=index_molecule[j]
            temp.append(index)
            MolSet[index].index_droplet=i
            MolSet[index].barcode=barcode[i]
            totalseqlen=totalseqlen+MolSet[index].length
        droplet_container.append(temp)
    return MolSet

def child_initialize(_MolSet,_reflist):
     global MolSet
     global reflist
     MolSet = _MolSet
     reflist=_reflist

# I think this is the main program that runs stuff
# Yep, it is
def haploid(Par,lib):
    """This is the main program that runs the simulations by calling everything else"""
    global MolSet    
    global reflist
    droplet_container=[]
    (reflist,reftitle)=input_seq(Par.fastaHap1)
    if Par.hap==2:
        (reflist2,reftitle2)=input_seq(Par.fastaHap2)
        reflist.extend(reflist2)
        reftitle.extend(reftitle2)
    
    #recode cut position of long fragment
    #print('read template finished (library '+lib+')')
    MolSet=randomlong(Par,reflist,reftitle)
    #print('generate molecule finished (library '+lib+')')
    #calculate number of droplet
    assign_drop=deternumdroplet(N_frag,Par.molPerDroplet)
    #print('assign molecule to droplet finished (library '+lib+')')
    MolSet=selectbarcode(Par.barcodePool,assign_drop,MolSet,droplet_container)
    #print('assign barcode to molecule finished (library '+lib+')')
    #print('begin to simulate short reads, please wait...')
    pool = multiprocessing.Pool(int(Par.threads),initializer= child_initialize,initargs = (MolSet,reflist,))
    Mol_process=[]
    maxprocessor=int(len(MolSet)/int(Par.threads))
    t=0
    while t<len(MolSet):
        Mol_process.append(t)
        t += maxprocessor
    Mol_process.append(len(MolSet)-1)
    # parallelize the process
    for m in range(len(Mol_process)-1):
        pool.apply_async(SIMSR,(Mol_process[m],Mol_process[m+1],Par,lib,m,))
    pool.close()
    pool.join()
    out_f = snakemake.output["fw"]
    out_r = snakemake.output["rv"]
    os.system(f'touch {out_f}')
    os.system(f'touch {out_r}')
    for m in range(len(Mol_process)-1):
        os.system(f'cat {lib}_S1_L001_id{m}_R1_001.fq.gz >> {out_f} && rm {lib}_S1_L001_id{m}_R1_001.fq.gz')
        os.system(f'cat {lib}_S1_L001_id{m}_R2_001.fq.gz >> {out_r} && rm {lib}_S1_L001_id{m}_R2_001.fq.gz')
    #print('Library '+lib+' simulation completed!')
    return

def reverseq(seq):
    complementary=''
    rev_complementary=''
    for i in range(len(seq)):
        if seq[i]=='A':
           complementary+='T'
        elif seq[i]=='T':
           complementary+='A'
        elif seq[i]=='C':
           complementary+='G'
        elif seq[i]=='G':
           complementary+='C'
        elif seq[i]=='N':
           complementary+='N'
    rev_complementary=complementary[::-1]
    return rev_complementary

def Input_BarcodeQual(Par):
    with open(Par.barcodeQualityFile,"r") as f:
        line_index=0
        position=0
        Qual_dict=defaultdict(list)
        Prob_dict=defaultdict(list)
        for line in f:
            if line_index>0:
                linequal=line.strip('\t,\n')
                qualarray=linequal.split('\t')
                Qual_dict[qualarray[0]].append(ord(qualarray[1]))
                Prob_dict[qualarray[0]].append(float(qualarray[2]))
            line_index += 1
    return Qual_dict,Prob_dict

def Input_SeqQual(Par):
    with open(Par.seqQualityFile,"r") as f:
        line_index=0
        position=0
        Qual_dict=defaultdict(list)
        Prob_dict=defaultdict(list)
        Substitute_dict=defaultdict(list)
        for line in f:
            if line_index>0:
                change=[]
                linequal=line.strip('\t,\n')
                qualarray=linequal.split('\t')
                Qual_dict[qualarray[0]].append(ord(qualarray[1]))
                Prob_dict[qualarray[0]].append(float(qualarray[2]))
                Substitute_dict[(qualarray[0],ord(qualarray[1]))]=list(map(float,qualarray[3:]))
            line_index += 1
    return Qual_dict,Prob_dict,Substitute_dict

def SIMSR(start,end,Par,lib,jobid):
    MolSet_cand=MolSet[start:end]
    f_reads1 = gzip.open(lib+'_S1_L001_id'+str(jobid)+'_R1_001.fastq.gz',"wb")
    f_reads2 = gzip.open(lib+'_S1_L001_id'+str(jobid)+'_R2_001.fastq.gz',"wb")
    [SeqQual_dict,SeqProb_dict,SeqSubstitute_dict]=Input_SeqQual(Par)
    [BarcodeQual_dict,BarcodeProb_dict]=Input_BarcodeQual(Par)
    last_reads=0
    for i in range(len(MolSet_cand)):
        seq=reflist[MolSet_cand[i].seqidx]
        Seq_rand_qual=[]
        Barcode_rand_qual=[]
        num_reads=int(int(MolSet_cand[i].length/(Par.lenShortRead*2))*Par.coverageShortRead)
        if num_reads==0:
            continue
        All_forward=seq[MolSet_cand[i].start:MolSet_cand[i].end].upper()
        All_reverse=reverseq(All_forward)
        insert_size=np.random.normal(loc=Par.avgInsertShortRead, scale=Par.stdInsertShortRead,size=num_reads)
        Totalreads=[]
        new_reads=last_reads+num_reads
        Seq_new_qual=np.zeros((num_reads*2,Par.lenShortRead),dtype=int)
        Barcode_new_qual=np.zeros((num_reads,16),dtype=int)
        for m in range(16):
            Barcode_coll_phred=BarcodeQual_dict[str(m)]
            Barcode_coll_prob=BarcodeProb_dict[str(m)]
            Barcode_new_qual[:,m]=np.random.choice(Barcode_coll_phred,p=Barcode_coll_prob,size=(num_reads))
        for m in range(Par.lenShortRead):
            Seq_coll_phred=SeqQual_dict[str(m)]
            Seq_coll_prob=SeqProb_dict[str(m)]
            Seq_new_qual[:,m]=np.random.choice(Seq_coll_phred,p=Seq_coll_prob,size=(num_reads*2))
        Seq_new_qual1=Seq_new_qual[0:num_reads,:]
        Seq_new_qual2=Seq_new_qual[num_reads:2*num_reads,:]
        for j in range(num_reads):
            PE_read=pairend(Par,insert_size,MolSet_cand[i],Barcode_new_qual,Seq_new_qual1,Seq_new_qual2,All_forward,All_reverse,j,SeqSubstitute_dict)
            Totalreads.append(PE_read)
        for j in range(len(Totalreads)):
            read1seq=Totalreads[j].seq1
            read1qual=Totalreads[j].qual1
            read2seq=Totalreads[j].seq2
            read2qual=Totalreads[j].qual2
            read1N=read1seq[23:Par.lenShortRead]
            if read1N.count('N')>(Par.lenShortRead-23)*0.1 or read2seq.count('N')>Par.lenShortRead*0.1:
               continue
            readname='@ST-K00126:'+str(i+1)+':H5W53BBXX:'+str(MolSet_cand[i].start)+':'+str(MolSet_cand[i].end)+':'+str(Totalreads[j].start1)+':'+str(Totalreads[j].end1)
            f_reads1.write((readname+' 1:N:0\n').encode('utf-8'))
            f_reads1.write((read1seq+'\n').encode('utf-8'))
            f_reads1.write(('+\n').encode('utf-8'))
            f_reads1.write((read1qual+'\n').encode('utf-8'))
            #f_reads1.write(b'\n')
            f_reads2.write((readname+' 3:N:0\n').encode('utf-8'))
            f_reads2.write((read2seq+'\n').encode('utf-8'))
            f_reads2.write(('+\n').encode('utf-8'))
            f_reads2.write((read2qual+'\n').encode('utf-8'))
    f_reads1.close()
    f_reads2.close()
    return None

def pairend(Par,insert_size,MolSetX,Barcode_rand_qual,Seq_rand_qual1,Seq_rand_qual2,All_forward,All_reverse,index,SeqSubstitute_dict):
    is_read=int(np.absolute(insert_size[index]))
    start_for=int(np.random.uniform(low=1,high=MolSetX.length-Par.lenShortRead-1))
    if start_for+is_read+1>MolSetX.length:
       is_read=MolSetX.length-start_for-1
    end_for=start_for+int(is_read)
    forward_seq=All_forward[start_for:end_for]
    start_rev=MolSetX.length-end_for
    end_rev=start_rev+int(is_read)
    reverse_seq=All_reverse[start_rev:end_rev]
    read1=forward_seq[23:Par.lenShortRead]
    read2=reverse_seq[0:Par.lenShortRead]
    read1seq=''
    read2seq=''
    read1qual=''
    read2qual=''
    readerror1=np.random.choice([0,1],p=[1-Par.errorRate,Par.errorRate],size=(Par.lenShortRead-23))
    readerror2=np.random.choice([0,1],p=[1-Par.errorRate,Par.errorRate],size=(Par.lenShortRead))
    error1=readerror1[0:Par.lenShortRead].nonzero()
    error2=readerror2[0:Par.lenShortRead].nonzero()
    read1new=list(read1)
    read2new=list(read2)
    for i in error1[0]:
        if read1[i]=='A':
            rand_nuc=np.random.choice(['C','G','T','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual1[index,i])][0:4]))
        elif read1[i]=='C':
            rand_nuc=np.random.choice(['A','G','T','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual1[index,i])][4:8]))
        elif read1[i]=='G':
            rand_nuc=np.random.choice(['A','C','T','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual1[index,i])][8:12]))
        elif read1[i]=='T':
            rand_nuc=np.random.choice(['A','C','G','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual1[index,i])][12:16]))
        elif read1[i]=='N':
            rand_nuc=np.random.choice(['A','C','G','T'],1,p=np.asarray([0.25,0.25,0.25,0.25]))
        read1new[i]=rand_nuc[0]

    for i in error2[0]:
        if read2[i]=='A':
            rand_nuc=np.random.choice(['C','G','T','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual2[index,i])][0:4]))
        elif read2[i]=='C':
            rand_nuc=np.random.choice(['A','G','T','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual2[index,i])][4:8]))
        elif read2[i]=='G':
            rand_nuc=np.random.choice(['A','C','T','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual2[index,i])][8:12]))
        elif read2[i]=='T':
            rand_nuc=np.random.choice(['A','C','G','N'],1,p=np.asarray(SeqSubstitute_dict[(str(i),Seq_rand_qual2[index,i])][12:16]))
        elif read2[i]=='N':
            rand_nuc=np.random.choice(['A','C','G','T'],1,p=np.asarray([0.25,0.25,0.25,0.25]))
        read2new[i]=rand_nuc[0]
        #TODO REMOVE THE NNNNNNNN
    read1seq=MolSetX.barcode+'NNNNNNN'+''.join(read1new)
    read2seq=''.join(read2new)
    read1qual=''.join(map(chr,Barcode_rand_qual[index,:]))+'KKKKKKK'+''.join(map(chr,Seq_rand_qual1[index,23:Par.lenShortRead]))
    read2qual=''.join(map(chr,Seq_rand_qual2[index,:]))
    return Short_reads_PE(read1seq,read1qual,start_for,end_for,read2seq,read2qual,start_rev,end_rev)


Par = parameter()
Par.fastaHap1 = snakemake.input["genome1"]
Par.fastaHap1 = snakemake.input["genome2"]
Par.hap = snakemake.params["haplotypes"]
Par.avgInsertShortRead = snakemake.params["short_insert"]
Par.avgLenLongFrag = snakemake.params["long_len"]
Par.barcodePool = snakemake.input["barcodes"]
Par.barcodeQualityFile = snakemake.input["bc_quality"]
Par.coverageLongFrag = snakemake.params["long_coverage"]
Par.coverageShortRead = snakemake.params["short_coverage"]
Par.errorRate = snakemake.params["error_rate"]
Par.lenShortRead = snakemake.params["short_len"]
Par.molPerDroplet = snakemake.params["mol_per_droplet"]
Par.seqQualityFile = snakemake.input["seq_quality"]
Par.stdInsertShortRead = snakemake.params["short_insert_sd"]
Par.threads = snakemake.threads
if Par.hap == 2:
    Par.coverageLongFrag /= 2
    #TODO output prefix?
haploid(Par, snakemake.params["tempdir"])

