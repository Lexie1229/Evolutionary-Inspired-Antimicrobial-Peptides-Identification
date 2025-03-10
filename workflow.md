# Evolutionary-Inspired-Antimicrobial-Peptides-Identification

## Identification Workflow

### Paenibacillus Genomes Collection

* Collect detailed information about Paenibacillus genomes

```bash
mkdir -p ~/biodata/new_polymyxin/genome
cd ~/biodata/new_polymyxin/genome
mkdir -p complete_taxon uncomplete_taxon complete_untaxon uncomplete_untaxon

# Update data to 23/8/8
nwr download
nwr txdb
nwr ardb
nwr ardb --genbank

# Species information for the genus Paenibacillus
nwr info "Paenibacillus"
nwr member 44249

# Genomes: Paenibacillus genus bacteria
echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN (44249)
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >info_all_genome.tsv

# Genomes: assembly level was complete genome and chromosome, and the taxonomic information was determined
echo "
.headers ON
    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN (44249)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >info_complete_taxon_genome.tsv

# Genomes: assembly level was scaffold and contig, and the taxonomic information was determined
echo "
.headers ON
    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN (44249)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Contig', 'Scaffold')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >info_uncomplete_taxon_genome.tsv

# Genomes: assembly level was complete genome and chromosome, and the taxonomic information was uncertain
echo "
.headers ON
    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN (44249)
        AND species LIKE '% sp.%'
        AND organism_name LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >info_complete_untaxon_genome.tsv

# Genomes: assembly level was scaffold and contig, and the taxonomic information was uncertain
echo "
.headers ON
    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN (44249)
        AND species LIKE '% sp.%'
        AND organism_name LIKE '% sp.%'
        AND assembly_level IN ('Contig', 'Scaffold')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >info_uncomplete_untaxon_genome.tsv
```

* Collect download information about Paenibacillus genomes

```bash
# Genomes: assembly level was complete genome and chromosome, and the taxonomic information was determined
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND genus_id IN (44249)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite >download_complete_taxon_genome.tsv

# Genomes: assembly level was scaffold and contig, and the taxonomic information was determined
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND genus_id IN (44249)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND name NOT LIKE '% sp.%'
        AND assembly_level IN ('Contig', 'Scaffold')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite >download_uncomplete_taxon_genome.tsv


# Genomes: assembly level was complete genome and chromosome, and the taxonomic information was uncertain
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND genus_id IN (44249)
        AND species LIKE '% sp.%'
        AND organism_name LIKE '% sp.%'
        AND name LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite >download_complete_untaxon_genome.tsv

# Genomes: assembly level was scaffold and contig, and the taxonomic information was uncertain
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND genus_id IN (44249)
        AND species LIKE '% sp.%'
        AND organism_name LIKE '% sp.%'
        AND name LIKE '% sp.%'
        AND assembly_level IN ('Contig', 'Scaffold')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite >download_uncomplete_untaxon_genome.tsv
```

* Modify the name of Paenibacillus genus strains

```bash
cd ~/biodata/new_polymyxin/genome

# Abbreviate strain name
for name in complete_taxon uncomplete_taxon complete_untaxon uncomplete_untaxon; do
    cat "download_${name}_genome.tsv" |
        grep -v '^#' |
        tsv-uniq |
        perl ~/biodata/new_polymyxin/script/abbr_name.pl -c "1,2,3" -s '\t' -min 5 --shortsub |
        (echo -e '#name\tftp_path\torganism\tassembly_level' && cat ) |
        perl -nl -a -F"," -e '
            BEGIN{my %seen};
            /^#/ and print and next;
            /^organism_name/i and next;
            $seen{$F[3]}++; # ftp_path
            $seen{$F[3]} > 1 and next;
            $seen{$F[5]}++; # abbr_name
            $seen{$F[5]} > 1 and next;
            printf qq{%s\t%s\t%s\t%s\n}, $F[5], $F[3], $F[1], $F[4];
            ' |
        keep-header -- sort -k3,3 -k1,1 >abbr_${name}_genome.tsv
done
```

* Download Paenibacillus genomes

```bash
cd ~/biodata/new_polymyxin/genome

# Prepare
for name in complete_taxon uncomplete_taxon complete_untaxon uncomplete_untaxon; do
    nwr assembly abbr_${name}_genome.tsv -o $name
done

# Download
for name in complete_taxon uncomplete_taxon complete_untaxon uncomplete_untaxon; do
    bash $name/rsync.sh
done

# Check
for name in complete_taxon uncomplete_taxon complete_untaxon uncomplete_untaxon; do
    bash $name/check.sh
    rm $name/check.list
done
```

### Represent Sequences and pHMMs of Domains

* Representative sequences obtained by clustering

```bash
cd ~/biodata/new_polymyxin/search
mkdir -p 2refer/1domain/1raw 2refer/1domain/2cdhit

# Sequence of domains
cp ../phylogenetic/1domain/2record/*.fa 2refer/1domain/1raw/

# Clustering based on identity threshold
cd-hit -i 2refer/1domain/1raw/adomain.fa -o 2refer/1domain/2cdhit/adomain.fa -c 0.95 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/cdomain.fa -o 2refer/1domain/2cdhit/cdomain.fa -c 0.96 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/tdomain.fa -o 2refer/1domain/2cdhit/tdomain.fa -c 0.94 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/edomain.fa -o 2refer/1domain/2cdhit/edomain.fa -c 0.93 -n 5 -M 10000 -T 8 -d 0
cd-hit -i 2refer/1domain/1raw/tedomain.fa -o 2refer/1domain/2cdhit/tedomain.fa -c 0.96 -n 5 -M 10000 -T 8 -d 0
```

* Build pHMM models

```bash
cd ~/biodata/new_polymyxin/search
mkdir -p 2refer/2hmm/1raw

# Sequence of domains
cp ../phylogenetic/1domain/2record/*.fa 2refer/2hmm/1raw/

# Sequences deduplicate
for domain in adomain cdomain tdomain edomain tedomain; do
    mkdir -p 2refer/2hmm/2class/${domain}

    # Obtain domain type list
    cat 2refer/2hmm/1raw/${domain}.fa |
        grep '>' |
        sed 's/>//g' |
        cut -d'.' -f 1,2,3 |
        sort | uniq >2refer/2hmm/2class/${domain}.tsv

    # Classify domain sequences
    cat 2refer/2hmm/2class/${domain}.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 8 "
            cat 2refer/2hmm/1raw/${domain}.fa |
                grep -A 1 {} |
                sed 's/^--$//g' |
                grep -v '^$'  >2refer/2hmm/2class/${domain}/{}.fa
        "

    # Sequences deduplicate
    cat 2refer/2hmm/2class/${domain}.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 8 "
            cd-hit -i 2refer/2hmm/2class/${domain}/{}.fa -o 2refer/2hmm/2class/${domain}/{}.de.fa -c 1.00 -n 5 -M 10000 -T 8 -d 0
        "
done

# Build profile HMMs
for domain in adomain cdomain tdomain edomain tedomain; do
    mkdir -p 2refer/2hmm/3phmm/${domain}
    cat 2refer/2hmm/2class/${domain}.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 8 "
            mafft --auto 2refer/2hmm/2class/${domain}/{}.de.fa >2refer/2hmm/3phmm/${domain}/{}.mafft.fa
            trimal -in 2refer/2hmm/3phmm/${domain}/{}.mafft.fa -out 2refer/2hmm/3phmm/${domain}/{}.trim.fa
            esl-reformat stockholm 2refer/2hmm/3phmm/${domain}/{}.trim.fa >2refer/2hmm/3phmm/${domain}/{}.sto
            hmmbuild 2refer/2hmm/3phmm/${domain}/{}.hmm 2refer/2hmm/3phmm/${domain}/{}.sto
    "
done

# Build pHMM database
cat 2refer/2hmm/3phmm/*/*.hmm >2refer/all.hmm
hmmpress 2refer/all.hmm
```

### Representative Sequences Search 

* Use miniprot to search

```bash
cd ~/biodata/new_polymyxin/search

# miniprot index and search
mkdir -p 3search/1miniprot/1index 3search/1miniprot/2search
cat strains.all.tsv | parallel --no-run-if-empty --linebuffer -k -j 8 '
    miniprot -t 2 -d 3search/1miniprot/1index/{}.mpi 1genome/{}.fna
    miniprot -t 2 --outs=0.85 --gff --aln --trans 3search/1miniprot/1index/{}.mpi 2refer/represent_domain.all.fa >3search/1miniprot/2search/{}.gff3
'

# Extract miniprot search results
mkdir -p 3search/1miniprot/3result
cat strains.all.tsv | parallel --no-run-if-empty --linebuffer -k -j 8 '
    cat 3search/1miniprot/2search/{}.gff3 | 
        grep "##PAF" | 
        cut -f 2,3,6,7,8,9,10,15 |
        sed "s/ms:i://g" |
        sort -k4,4 -k6,6n -k7,7n -k8,8nr |
        uniq >3search/1miniprot/3result/{}.tsv
'

# Handle miniprot results: remove outliers and merge identical matches
cat strains.all.tsv | parallel --no-run-if-empty --linebuffer -k -j 8 '
    python 0script/remove_outliers.py -m miniprot -i 3search/1miniprot/3result/{}.tsv -o 3search/1miniprot/4filter/{}.outliers.tsv
    python 0script/combine_match.py -t 0.95 -i 3search/1miniprot/4filter/{}.outliers.tsv -o 3search/1miniprot/5site/{}.tsv
'
```

* Use blast to search

```bash
cd ~/biodata/new_polymyxin/search

# blast creatdb and search
mkdir -p 3search/2blast/2search
cat strains.all.tsv | parallel --no-run-if-empty --linebuffer -k -j 2 '
    mkdir -p 3search/2blast/1db/{}
    makeblastdb -in 1genome/{}.fna -dbtype nucl -out 3search/2blast/1db/{}/{} 1>/dev/null
    tblastn -query 2refer/represent_domain.all.fa -db 3search/2blast/1db/{}/{} -out 3search/2blast/2search/{}.tsv -num_threads 4 -outfmt "6 qseqid qlen sstrand sseqid slen sstart send bitscore pident evalue qcovs"
    sed -i "s/plus/+/g;s/minus/-/g" 3search/2blast/2search/{}.tsv
'

# Extract blast search results and Format results
mkdir -p 3search/2blast/3result
cat strains.all.tsv | parallel --no-run-if-empty --linebuffer -k -j 8 '
    perl 0script/blast2miniprot.pl 3search/2blast/2search/{}.tsv >3search/2blast/3result/{}.tsv
'

# Handle blast results: remove outliers, matches with low confidence, matches with low quality, and merge identical matches
mkdir -p 3search/2blast/4filter 3search/2blast/5site
cat strains.all.tsv | parallel --no-run-if-empty --linebuffer -k -j 8 '
    python 0script/remove_outliers.py -m blast -i 3search/2blast/3result/{}.tsv -o 3search/2blast/4filter/{}.outliers.tsv
    cat 3search/2blast/4filter/{}.outliers.tsv |
        tsv-filter --ge 9:70 --le 10:1.0e-20 --ge 11:50 |
        sort -k4,4 -k6,6n -k7,7n -k8,8nr |
        uniq >3search/2blast/4filter/{}.confidence.tsv
    python 0script/quality_control.py 3search/2blast/4filter/{}.confidence.tsv 3search/2blast/4filter/{}.quality.tsv
    python 0script/combine_match.py -t 0.95 -i 3search/2blast/4filter/{}.quality.tsv -o 3search/2blast/5site/{}.tsv
'
```

### Extract Matched Domains

* Extracts sequences of matched domains

```bash
cd ~/biodata/new_polymyxin/search

# Merge search results of miniprot and blast
mkdir -p 3search/3site/1raw
cat strains.all.tsv | parallel --no-run-if-empty --linebuffer -k -j 8 '
    cat 3search/1miniprot/5site/{}.tsv 3search/2blast/5site/{}.tsv |
        awk -F"\t" "{ if (\$5 == 0) { \$5 = 1; } print \$0; }" OFS="\t" |
        sort -k3,3 -k5,5n -k6,6n |
        uniq > 3search/3site/1raw/{}.tsv
'

# Merge site information and divide the cluster
mkdir -p 3search/3site/2site 3search/3site/3cluster
cat strains.all.tsv | parallel --no-run-if-empty --linebuffer -k -j 8 '
    python 0script/recombine_site.py -t 0.95 -i 3search/3site/1raw/{}.tsv -o 3search/3site/2site/{}.tsv
    python 0script/divide_cluster.py 3search/3site/2site/{}.tsv 3search/3site/3cluster/{}.tsv
'

# Extract site information of matched Domains
mkdir -p 3search/4subseq/1region
cat strains.all.tsv | parallel --no-run-if-empty --linebuffer -k -j 8 '
    cat 3search/3site/3cluster/{}.tsv |
        grep -v ">" |
        cut -f 7,9,10 |
        perl -nla -e "
            my (\$contig, \$start, \$end) = split /\t/;
            print qq{\$contig:\$start-\$end};
        " >3search/4subseq/1region/{}.tsv
'

# Extract sequences of matched Domains
while IFS= read -r strain; do
    mkdir -p 3search/4subseq/2nucl/${strain}
    cat 3search/4subseq/1region/${strain}.tsv | parallel --colsep '\n' --no-run-if-empty --linebuffer -k -j 8 "
        faops region -l 0 1genome/${strain}.fna <(echo {}) 3search/4subseq/2nucl/${strain}/{}.fa
"
done <strains.all.tsv

# Six frame translation
while IFS= read -r strain; do
    mkdir -p 3search/4subseq/3pro/${strain}
    cat 3search/4subseq/1region/${strain}.tsv | parallel --colsep '\n' --no-run-if-empty --linebuffer -k -j 8 "
        perl 0script/six_frame_translation.pl 3search/4subseq/2nucl/${strain}/{}.fa 3search/4subseq/3pro/${strain}/{}.fa
        "
done <strains.all.tsv
```

### Determination of Domain Types

```bash
cd ~/biodata/new_polymyxin/search

# hmmscan based on pHMMs
while IFS= read -r strain; do
    mkdir -p 3search/5hmmscan/1scan/${strain}
    cat 3search/4subseq/1region/${strain}.tsv | parallel --colsep '\n' --no-run-if-empty --linebuffer -k -j 2 "
            hmmscan --cpu 4 --tblout 3search/5hmmscan/1scan/${strain}/{}.txt 2refer/all.hmm 3search/4subseq/3pro/${strain}/{}.fa
        "
done <strains.all.tsv

# Extract hmmscan results
while IFS= read -r strain; do
    mkdir -p 3search/5hmmscan/2result/${strain}
    cat 3search/4subseq/1region/${strain}.tsv | parallel --colsep '\n' --no-run-if-empty --linebuffer -k -j 8 "
            python 0script/extract_hmm.py 3search/5hmmscan/1scan/${strain}/{}.txt 3search/5hmmscan/2result/${strain}/{}.tsv
        "
done <strains.all.tsv

# Merge different results reflecting the same site
while IFS= read -r strain; do
    mkdir -p 3search/5hmmscan/3merge/${strain}
    cat 3search/4subseq/1region/${strain}.tsv | parallel --colsep '\n' --no-run-if-empty --linebuffer -k -j 8 "
        if [ -s 3search/5hmmscan/2result/${strain}/{}.tsv ]; then
            python 0script/merge_lines.py -i 3search/5hmmscan/2result/${strain}/{}.tsv -o 3search/5hmmscan/3merge/${strain}/{}.tsv
        fi
    "
done <strains.all.tsv

# Merge hmmscan results
mkdir -p 3search/5hmmscan/4combine
cat strains.all.tsv | parallel --colsep '\n' --no-run-if-empty --linebuffer -k -j 8 '
    if [ "$(ls -A 3search/5hmmscan/3merge/{}/ 2>/dev/null)" ]; then
        cat 3search/5hmmscan/3merge/{}/*.tsv >3search/5hmmscan/4combine/{}.tsv
    else
        touch 3search/5hmmscan/4combine/{}.tsv
    fi
'

# Merge search and hmmscan results
mkdir -p 3search/5hmmscan/5all
cat strains.all.tsv | parallel --colsep '\n' --no-run-if-empty --linebuffer -k -j 8 '
    python 0script/combine_result.py -s 3search/3site/3cluster/{}.tsv -m 3search/5hmmscan/4combine/{}.tsv -o 3search/5hmmscan/5all/{}.tsv
'
```

## Build Protein Trees of Domains

* Build protein trees based on different domains

```bash
cd ~/biodata/new_polymyxin/phylogenetic

# Build protein trees
for domain in adomain cdomain tdomain edomain tedomain; do 
    mkdir -p 1domain/3tree/${domain}
    # Multiple sequence alignment
    mafft --auto --thread 4 1domain/2record/${domain}.fa >1domain/3tree/${domain}/${domain}.mafft.fa
    # Trim
    trimal -in 1domain/3tree/${domain}/${domain}.mafft.fa -out 1domain/3tree/${domain}/${domain}.trim.fa
    # Build trees
    FastTree 1domain/3tree/${domain}/${domain}.trim.fa >1domain/3tree/${domain}/${domain}.newick
done

# Sort nodes of protein trees
for domain in adomain cdomain tdomain edomain tedomain; do 
    nwr order --nd 1domain/3tree/${domain}/${domain}.newick -o 1domain/3tree/${domain}/${domain}.order.newick
done
```

## Build Species Tree

* Build species tree based on 120 single-copy proteins(bac120)

```bash
cd ~/biodata/new_polymyxin/phylogenetic/5species

# download bac120 related list and profile HMMs
nwr kb bac120 -o HMM
cp HMM/bac120.lst HMM/marker.lst

# Extract genome protein file - paenibacillus
for strain in $(cat strain.tsv); do 
    mkdir -p Protein/${strain}
    find ../../genome/*/${strain} -name "*_protein.faa.gz" | 
        parallel --no-run-if-empty --linebuffer -k -j 8 "
            cp {} Protein/${strain}/pro_faa.fa.gz; 
        "
done

# Extract genome protein file - outgroup
for strain in $(cat outgroup.tsv); do
    mkdir -p Protein/${strain}
    cp outgroup/${strain}.faa Protein/${strain}/pro_faa.fa; 
    gzip Protein/${strain}/pro_faa.fa
done

# combine protein file
cat strain.tsv outgroup.tsv >all.strain.tsv

# Extract bac120 list
cat all.strain.tsv |
    while read STRAINS; do
        if [[ -s Protein/"${STRAINS}"/bac120.tsv ]]; then
            continue
        fi
        if [[ ! -f Protein/"${STRAINS}"/pro_faa.fa.gz ]]; then
            continue
        fi

        echo >&2 "${STRAINS}"

        cat HMM/marker.lst |
            parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 "
                gzip -dcf Protein/${STRAINS}/pro_faa.fa.gz |
                    hmmsearch --cut_nc --noali --notextw HMM/hmm/{}.HMM - |
                    grep '>>' |
                    perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), q({}), \$1; '
            " >Protein/${STRAINS}/bac120.tsv
    done

# Formate bac120 name
mkdir -p Domain
cat all.strain.tsv |
    while read STRAINS; do
        cat Protein/${STRAINS}/bac120.tsv |
            env STRAINS="${STRAINS}" perl -F'\t' -lane '
                $bac120 = $F[0];
                $wp = $F[1];
                $name = "$ENV{STRAINS}+$bac120+$wp";
                print "$name\t$ENV{STRAINS}\t$bac120\t$wp";
            ' 
    done >Domain/bac120.tsv

# Extract bac120 sequences
cat all.strain.tsv |
    while read STRAINS; do
        cat Protein/${STRAINS}/bac120.tsv | 
            parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 "
                hnsm some Protein/${STRAINS}/pro_faa.fa.gz <(echo {2}) | sed 's/>.*$/>${STRAINS}+{1}+{2}/g'
        "
    done |
        hnsm dedup stdin |
        hnsm gz stdin -o Domain/bac120.fa

# Extract bac120 sequence for each genome
cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        echo >&2 "==> marker [{}]"
            
        mkdir -p Domain/{}

        hnsm some Domain/bac120.fa.gz <(
            cat Domain/bac120.tsv |
                tsv-filter --str-eq "3:{}" |
                tsv-select -f 1 |
                rgr dedup stdin
            ) >Domain/{}/{}.pro.fa
    '

# Multiple sequence alignment
cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        echo >&2 "==> marker [{}]"

        if [ ! -s Domain/{}/{}.pro.fa ]; then
            exit
        fi
        if [ -s Domain/{}/{}.aln.fa ]; then
            exit
        fi

        mafft --auto Domain/{}/{}.pro.fa >Domain/{}/{}.aln.fa
    '

# Tandem bac120 protein
cat HMM/marker.lst |
    while read marker; do
        if [ ! -s Domain/${marker}/${marker}.pro.fa ]; then
            continue
        fi
        if [ ! -s Domain/${marker}/${marker}.aln.fa ]; then
            continue
        fi

        cat Domain/${marker}/${marker}.aln.fa | faops filter -l 0 stdin stdout

        # empty line for .fas
        echo
    done >Domain/bac120.aln.fas
sed -i 's/+/ /g' Domain/bac120.aln.fas
cat Domain/bac120.tsv |
    cut -f 2 |
    rgr dedup stdin |
    sort |
    fasops concat Domain/bac120.aln.fas stdin -o Domain/bac120.aln.fa

# Trim sequences
trimal -in Domain/bac120.aln.fa -out Domain/bac120.trim.fa -automated1

# Build species trees
FastTree -fastest -noml Domain/bac120.trim.fa >Domain/bac120.trim.newick

# Definite root
mkdir -p tree
nw_reroot Domain/bac120.trim.newick Brevibac_agri_GCF_004117055_1 Brevibac_brevis_GCF_900637055_1 |
    nwr order stdin --nd --an >tree/bac120.reroot.newick
```

## Protein Structure Prediction and Comparison

* Prepare json file for Alphafold

```bash
cd ~/biodata/new_polymyxin/alphafold

# Format conversion: fasta → json
for domain in cdomain tdomain; do 
    for file in $(find ${domain} -type f -name "*.fa"); do
        name=$(basename ${file} .fa) 
        python alphafold_fasta2json.py -i ${domain}/${name}.fa -o ${domain}/${name}.json
    done
done
```

* Structure comparison: based on TM-score

```bash
cd ~/biodata/new_polymyxin/alphafold

# Calculate TM-score
for domain in cdomain tdomain; do
    # Extract cif file: model 0
    python alphafold_cif.py -i ${domain}_result -o ${domain}_cif

    # Calculate TM-score
    USalign -infmt1 3 -infmt2 3 -outfmt 2 -dir ${domain}_cif/ ${domain}.result.tsv >${domain}.score.tsv

    # Extract anno information
    python extract_anno.py -i ${domain}.name.tsv -o ${domain}.anno.tsv

    # Calculate TM-score mean
    python extract_score_mean.py -i ${domain}.score.tsv -o ${domain}.score.mean.tsv
done

# Result conversion: TM-score → distance
for domain in cdomain tdomain; do
    sed '1d' ${domain}.score.mean.tsv | cut -f 1,2,6 >${domain}.tmp.tsv
    hnsm convert --mode matrix ${domain}.tmp.tsv >${domain}.distance.mean.tsv
    rm ${domain}.tmp.tsv
done

# Heatmap
Rscript heatmap_mean.R cdomain
Rscript heatmap_mean.R tdomain
```

