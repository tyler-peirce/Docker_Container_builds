using FASTX
using BioSequences

struct Exon
    id::String
    interval::UnitRange{Int}
end

struct GeneModel
    gene::String
    strand::Char
    exons::Vector{Exon}
end

function splice(genome::LongDNA, gm::GeneModel)
    sort!(gm.exons, by=x->x.id)
    spliced = LongDNA{4}()
    for e in gm.exons
        if e.interval.start < 0
            append!(spliced, genome[mod1(e.interval.start, length(genome)):length(genome)])
        end
        append!(spliced, genome[max(1, e.interval.start):min(length(genome), e.interval.stop)])
        if e.interval.stop > length(genome)
            append!(spliced, genome[1:mod1(e.interval.stop, length(genome))])
        end
    end

    #remove stop codon or partial stop codon
    spliced = spliced[1:end - mod1(length(spliced), 3)]

    return spliced
end

function get_attribute(key::AbstractString, attributestring::AbstractString)
    attributes = split(attributestring, ";")
    attribute = attributes[findfirst(x->startswith(x, key), attributes)]
    attribute[length(key)+2:end]
end

function rc(glength::Int, pos::Int)
    return glength - pos + 1
end
    
function main(fastadir::String, gffdir::String)
    fasta_files = sort!(filter!(x->endswith(x, ".fasta") || endswith(x, ".fa"), readdir(fastadir; join=true)))
    gff_files = sort!(filter!(x->endswith(x, ".gff3") || endswith(x, ".gff"), readdir(gffdir; join=true)))
    println("FASTA files found: ", fasta_files)
    println("GFF files found: ", gff_files)
    for g in gff_files
        f = joinpath(fastadir, basename(g)[1:end-4])
        !isfile(f) && continue
        println("Processing FASTA: ", f)
        println("Processing GFF: ", g)
        println(f)
        println(g)
        species_name = first(split(basename(f), ".f"))
        fwd = FASTA.Reader(open(f)) do infile
            first(infile)
        end
        fwdseq = FASTA.sequence(LongDNA{4}, fwd)
        rev = FASTA.Record(FASTA.identifier(fwd), reverse_complement(fwdseq))
        revseq = FASTA.sequence(LongDNA{4}, rev)
        gff_lines = open(g) do infile
            readlines(infile)
        end
        filter!(x->!startswith(x,"#"),gff_lines)
        #= source_line = gff_lines[findfirst(occursin.("\tsource\t",gff_lines))]
        fields = split(source_line, '\t')
        genome_length = parse(Int, fields[5])
        println(genome_length) =#

        genome_length = length(fwdseq)

        filter!(x -> occursin("\trRNA\t", x) || occursin("\tCDS\t", x), gff_lines)

        genemodels = GeneModel[]
        current_strand = '.'
        current_genemodel = ""
        current_exons = Exon[]
        for line in gff_lines
            fields = split(line, '\t')
            strand = fields[7][1]
            if get_attribute("Name", fields[9]) ≠ current_genemodel
                if current_genemodel ≠ ""
                    push!(genemodels, GeneModel(current_genemodel, current_strand, current_exons))
                end
                current_genemodel = get_attribute("Name", fields[9])
                current_strand = strand
                current_exons = Exon[]
            end
            if strand == '+'
                push!(current_exons, Exon(get_attribute("Name", fields[9]), parse(Int, fields[4]):parse(Int, fields[5])))
            else
                push!(current_exons, Exon(get_attribute("Name", fields[9]), rc(genome_length, parse(Int, fields[5])):rc(genome_length, parse(Int, fields[4]))))
            end
        end
        push!(genemodels, GeneModel(current_genemodel, current_strand, current_exons))
        println(pwd())
        if !isdir("cds")
            mkdir("cds")
        end
        if !isdir("proteins")
            mkdir("proteins")
        end
        for gm in genemodels
            startswith(gm.gene, "matK") && continue
            endswith(gm.gene, "-2") && continue
            if gm.strand == '+'
                cds = splice(fwdseq, gm)
            else
                cds = splice(revseq, gm)
            end
            FASTA.Writer(open("cds/" * gm.gene * ".fa", "a")) do outfile
                write(outfile, FASTA.Record(species_name, cds))
            end
            if mod(length(cds), 3) ≠ 0
                println("misannotated CDS: ", gm.gene)
            end
            protein = translate(cds, code=BioSequences.vertebrate_mitochondrial_genetic_code)
            FASTA.Writer(open("proteins/" * gm.gene * ".fa", "a")) do outfile
                write(outfile, FASTA.Record(species_name, protein))
            end
        end
    end
end

main(ARGS[1], ARGS[2])