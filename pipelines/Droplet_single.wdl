workflow main {
  String root
  String fastq1
  String fastq2
  String outdir
  String refdir
  String Rscript
  String Python3
  String ID
  String gtf
  String ?runID
  String config
  String ?lib
  Int ?expectCell
  Int ?forceCell
  Int ?umilow
  String ?species
  String ?original
  String ?SampleTime
  String ?ExperimentalTime
  call makedir {
    input:
    Dir=outdir,
  }
  call parseFastq {
    input:
    lib=lib,
    config=config,
    fastq1=fastq1,
    fastq2=fastq2,
    outdir=makedir.Outdir,
    runID=runID,
    root=root,
  }
  call fastq2bam {
    input:
    lib=lib,
    fastq=parseFastq.fastq,
    outdir=outdir,
    refdir=refdir,
    root=root
  }
  call sortBam {
    input:
    lib=lib,
    bam=fastq2bam.bam,
    gtf=gtf,
    root=root,
    outdir=outdir
  }
  call cellCount {
    input:
    lib=lib,
    bam=sortBam.anno,
    outdir=outdir,
    root=root,
    rawlist=parseFastq.rawlist
  }
  call cellCalling {
    input:
    lib=lib,
    root=root,
    count=cellCount.count,
    outdir=outdir,
    Rscript=Rscript,
    expectCell=expectCell,
    forceCell=forceCell,
    umilow=umilow
  }
  call countMatrix {
    input:
    lib=lib,
    root=root,
    list=cellCalling.list,
    outdir=outdir,
    anno=sortBam.anno,
    Rscript=Rscript,
    expectCell = expectCell,
    Python3=Python3,
    ID=ID
  }
  call report {
    input:
    ID=ID,
    root=root,
    outdir=outdir,
    species=species,
    Python3=Python3,
    original=original,
    SampleTime=SampleTime,
    ExperimentalTime=ExperimentalTime,
    matrix=countMatrix.matrix
  }
}

task makedir {
  String Dir
  command {
    if [ -f ${Dir}/symbol/makedir_sigh.txt ];then
      echo "makedir node success"
    else
      mkdir -p ${Dir}
      mkdir -p ${Dir}/outs
      mkdir -p ${Dir}/temp
      mkdir -p ${Dir}/report
      mkdir -p ${Dir}/symbol
      echo "[`date +%F` `date +%T`] workflow start" > ${Dir}/workflowtime.log
      echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${Dir}/symbol/makedir_sigh.txt
    fi
  }
  output {
    String Outdir="${Dir}"
  }
}

task parseFastq {
  String config
  String fastq1
  String fastq2
  String outdir
  String ?runID
  String root
  String ?lib
  Int cpp=1
  command {
    if [ -f ${outdir}/symbol/parseFastq_sigh.txt ];then
      echo "parseFastq node success"
    else
      if [ -f ${default=abjdbashj lib} ]; then
      source ${lib}
      fi
      ${root}/bin/PISA parse -t ${cpp} -f -q 20 -dropN -config ${config} -cbdis ${outdir}/temp/barcode_counts_raw.txt -run ${default=1 runID} -report ${outdir}/report/sequencing_report.csv ${fastq1} ${fastq2} |gzip -c  > ${outdir}/temp/reads.fq.gz
      head -n 50000 ${outdir}/temp/barcode_counts_raw.txt |cut -f1 > ${outdir}/temp/barcode_raw_list.txt
      echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/parseFastq_sigh.txt
    fi
  }
  output {
    String count="${outdir}/temp/barcode_counts_raw.txt"
    String fastq="${outdir}/temp/reads.fq.gz"
    String sequencingReport="${outdir}/report/sequencing_report.csv"
    String rawlist="${outdir}/temp/barcode_raw_list.txt"
  }
}

task fastq2bam {
  String fastq
  String outdir
  String refdir
  String root
  String ?lib
  Int cpp=4
  command {
    if [ -f ${outdir}/symbol/fastq2bam_sigh.txt ];then
      echo "fastq2bam node success"
    else
      if [ -f ${default=abjdbashj lib} ]; then
        source ${lib}
      fi
      ${root}/bin/STAR --outStd SAM --outSAMunmapped Within --runThreadN ${cpp} --genomeDir ${refdir} --readFilesCommand gzip -dc --readFilesIn ${fastq} --outFileNamePrefix ${outdir}/temp/ 1> aln.sam
      ${root}/bin/PISA sam2bam -@ ${cpp} -k -o ${outdir}/temp/aln.bam -report ${outdir}/report/alignment_report.csv aln.sam && rm -f aln.sam 
      echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/fastq2bam_sigh.txt
    fi
  }
  output {
    String bam="${outdir}/temp/aln.bam"
    String alnReport="${outdir}/report/alignment_report.csv"
  }
}

task sortBam {
  String bam
  String root
  String outdir
  String gtf
  String ?lib
  Int cpp=1
  command {
    if [ -f ${outdir}/symbol/sortBam_sigh.txt ];then
      echo "sortBam node success"
    else
      if [ -f ${default=abjdbashj lib} ]; then
        source ${lib}
      fi
      ${root}/bin/sambamba sort -t ${cpp} -o ${outdir}/temp/sorted.bam ${outdir}/temp/aln.bam
      ${root}/bin/PISA anno -gtf ${gtf} -o ${outdir}/temp/annotated.bam -report ${outdir}/report/annotated_report.csv ${outdir}/temp/sorted.bam
      ${root}/bin/PISA corr -tag UB -tags-block CB,GN -@ ${cpp} -t ${cpp} -o ${outdir}/outs/final.bam ${outdir}/temp/annotated.bam
      echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/sortBam_sigh.txt
    fi
  }
  output {
    String anno="${outdir}/temp/annotated.bam"
  }
}

task cellCount {
  String bam
  String outdir
  String root
  String rawlist
  String ?lib
  Int cpp=1
  command {
    if [ -f ${outdir}/symbol/cellCount_sigh.txt ];then
      echo "cellCount node success"
    else
      if [ -f ${default=abjdbashj lib} ]; then
        source ${lib}
      fi
      ${root}/bin/PISA attrcnt -cb CB -tags UB,GN -@ ${cpp} -dedup -o ${outdir}/temp/cell_stat.txt -list ${rawlist} ${bam}
      echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/cellCount_sigh.txt
    fi
  }
  output {
    String count="${outdir}/temp/cell_stat.txt"
  }
}
task cellCalling {
  String count
  String outdir
  String Rscript
  String root
  Int ?expectCell
  Int ?forceCell
  Int ?umilow
  String ?lib
  Int cpp=1
  command {
    if [ -f ${outdir}/symbol/cellCalling_sigh.txt ];then
      echo "cellCalling node success"
    else
      if [ -f ${default=abjdbashj lib} ]; then
        source ${lib}
      fi
      ${Rscript} ${root}/scripts/scRNA_cell_calling.R -i ${count} -o ${outdir}/outs -e ${default=0 expectCell}  -f ${default=0 forceCell} -l ${default=1000 umilow}
      echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/cellCalling_sigh.txt
    fi
  }
  output {
    String list="${outdir}/outs/cell_barcodes.txt"
  }
}
task countMatrix {
  String root
  String list
  String outdir
  String anno
  String ?lib
  String Rscript
  String ?Python3
  String ID
  Int ?expectCell
  Int cpp=1
  command {
    if [ -f ${outdir}/symbol/countMatrix_sigh.txt ];then
      echo "countMatrix node success"
    else
      if [ -f ${default=abjdbashj lib} ]; then
        source ${lib}
      fi
      ${root}/bin/PISA count -@ ${cpp} -tag CB -anno_tag GN -umi UB -o ${outdir}/outs/count_mtx.tsv -list ${list} ${anno}
      gzip ${outdir}/outs/count_mtx.tsv
      ${Python3} ${root}/scripts/scRNA_scanpy_clustering.py ${outdir}/outs/count_mtx.tsv.gz ${outdir}/report/
      echo "[`date +%F` `date +%T`] workflow end" >> ${outdir}/workflowtime.log
      echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/countMatrix_sigh.txt
    fi
  }
  output {
    String matrix = "${outdir}/outs/count_mtx.tsv.gz"
  }
}
task report {
  String Python3
  String ID
  String ?species
  String ?original
  String ?SampleTime
  String ?ExperimentalTime
  String outdir
  String root
  String matrix
  Int cpp=1
  command {
    if [ -f ${outdir}/symbol/report_sigh.txt ];then
      echo "report node success"
    else
      echo "Name,${ID}" > ${outdir}/report/sample.csv
      echo "Original,${default=Null original}" >>  ${outdir}/report/sample.csv
      echo "Species,${default=Null species}" >> ${outdir}/report/sample.csv
      echo "Sampling time,${default=Null SampleTime}" >> ${outdir}/report/sample.csv
      echo "Experimental time,${default=Null ExperimentalTime}" >> ${outdir}/report/sample.csv
      ${Python3} ${root}/scripts/idrop.py single ${outdir}/report iDrop_${ID}
      echo "rm ${outdir}/temp/*bam ${outdir}/temp/*fq" > ${outdir}/clean.sh
      echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/report_sigh.txt
    fi
  }
}

