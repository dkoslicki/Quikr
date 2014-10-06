# ==============================================================================
# ConvertToCAMIOutput.jl
#
# Authors: David Koslicki (david.koslicki@math.oregonstate.edu)
#
# Takes a raw reconstruction vector (flat text file of floats on the same basis
# as the given taxonomy file) and outputs the classification in the CAMI format.
# ==============================================================================

using ArgParse

#Parse arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input_file", "-i"
			help = "Input raw reconstruction text file"
        "--taxonomy_file", "-t"
			help = "taxonomy file, the ith line is the taxonomy of the ith training organism"
        "--output_level", "-l"
			help = "Output evolutionary related level, the lth entry of [1, .95, .8, .7, .6, .5, .4, .3, .2, .1]. Level of 0 means to sum over all of them."
			default = 0
        "--output_taxonomic_rank", "-r"
			help = "Output taxonomic rank, either an integer: 1, or a range [1:2]"
			default = "[1:7]"
		"--output_file", "-o"
			help = "Output text file"
		"--sample_ID", "-I"
			help = "Sample ID"
			default = "SAMPLEID"
		"--contestant_ID", "-C"
			help = "Contestant ID"
			default = "CONTESTANTID"
    end
    return parse_args(s)
end

#Parse the args
parsed_args = parse_commandline()
input_file = parsed_args["input_file"]
taxonomy_file = parsed_args["taxonomy_file"]
output_level = int(parsed_args["output_level"])
output_taxonomic_rank = eval(parse(parsed_args["output_taxonomic_rank"]))
output_file = parsed_args["output_file"]
sample_ID = parsed_args["sample_ID"]
contestant_ID = parsed_args["contestant_ID"]

#Read in the input file
fid = open(input_file,"r")
input = readlines(fid)
close(fid)
input = map(x->float(strip(x)), input)

#Next, read in the taxonomy file
fid = open(taxonomy_file,"r")
taxonomy = readlines(fid)
close(fid)
taxonomy = map(x->strip(split(x)[3]), taxonomy)
num_organisms = length(taxonomy)

#If the output_level is not 0, just do that level
if ~(output_level==0)

	#Just select the ones in the output level of interest
	indicies_of_interest = ((output_level-1) * num_organisms + 1):(output_level * num_organisms)

	#First, select the portion of the taxonomy that has a nonzero entry in the reconstruction
	cutoff = .00001
	support = indicies_of_interest[find(input[indicies_of_interest] .> cutoff)] #Support in the indicies of interest
	nonzero_taxonomy = taxonomy[(support .- (num_organisms*(output_level-1)))] #Shift everything left to the start since taxonomy has only num_organisms length
elseif output_level==0
	if ~(int(length(input)/num_organisms)==length(input)/num_organisms)
		error("Input reconstruction is not a multiple of the number of organisms in the taxonomy file")
	else
		#Then add up each chunk of size num_organisms
		input_temp = zeros(num_organisms);
		for output_level_temp = 1:int(length(input)/num_organisms)
			input_temp = input_temp + input[((output_level_temp-1) * num_organisms + 1):(output_level_temp * num_organisms)];
		end
		cutoff = .00001
		support = find(input_temp .> cutoff) #Support in the indicies of interest
		nonzero_taxonomy = taxonomy[support]
	end
else
	error("non-valid output_level")
end


#open the output file
output_file_handle = open(output_file,"w")

#Write the header
write(output_file_handle,"# CAMI Submission for Taxonomic Profiling\n")
write(output_file_handle, "@Task:Profiling\n")
write(output_file_handle,"@Version:1.0\n")
write(output_file_handle,"@ContestantID:CONTESTANTID\n")
write(output_file_handle,"@SampleID:$(sample_ID)\n")
write(output_file_handle,"@Ranks: superkingdom|phylum|class|order|family|genus|species|strain\n")
write(output_file_handle,"\n")
write(output_file_handle,"@@TAXID\tRANK\tTAXPATH\tTAXPATH_SN\tPERCENTAGE\n")

#Now for each of the ranks, get the unique names, then loop over the non-zero taxonomy, increasing the value of the unique taxa name, then output these to the file
if typeof(output_taxonomic_rank) == Int64
	taxa_rank_list = [output_taxonomic_rank]
elseif typeof(output_taxonomic_rank) == Array{Int64,1}
	taxa_rank_list = output_taxonomic_rank
else
	error("Input taxonomic rank should be an integer, or else don't include the option to output all ranks")
end

for taxa_rank = taxa_rank_list	
	taxa_names = cell(0);
	for taxonomy_string = nonzero_taxonomy
		nonzero_taxonomy_split = split(taxonomy_string,"|")
		if length(nonzero_taxonomy_split) >= taxa_rank
			taxa_name = join(nonzero_taxonomy_split[1:taxa_rank],"|")
			append!(taxa_names, {taxa_name})
		end
	end
	unique_taxa_names = sort(unique(taxa_names)); #This assumes that there's a bijection between taxa names and tax IDs

	#Now loop through each of the non_zero taxonomies, see if the taxa name matches, and then add this to the abundances
	taxa_abundances = Dict();
	for unique_taxa_name = unique_taxa_names
			taxa_abundances[unique_taxa_name] = 0
	end
	nonzero_taxonomy_counter = 1
	for taxonomy_string = nonzero_taxonomy
		nonzero_taxonomy_split = split(taxonomy_string,"|")
		if length(nonzero_taxonomy_split) >= taxa_rank
			taxa_name = join(nonzero_taxonomy_split[1:taxa_rank],"|")
			taxa_abundances[taxa_name] = taxa_abundances[taxa_name] + input[support[nonzero_taxonomy_counter]]
		end
		nonzero_taxonomy_counter = nonzero_taxonomy_counter + 1
	end
	for unique_taxa_name = unique_taxa_names
		taxID = split(split(unique_taxa_name,"|")[end],"_")[3];
		write(output_file_handle, "$(taxID)")
		write(output_file_handle, "\t")
		
		rankAbvr = split(split(unique_taxa_name,"|")[end],"_")[1];
		if rankAbvr == "k"
			rank = "superkingdom"
		elseif rankAbvr == "p"
			rank = "phylum"
		elseif rankAbvr == "c"
			rank = "class"
		elseif rankAbvr == "o"
			rank = "order"
		elseif rankAbvr == "f"
			rank = "family"
		elseif rankAbvr == "g"
			rank = "genus"
		elseif rankAbvr == "s"
			rank = "species"
		elseif rankAbvr == "t"
			rank = "string"
		else
			rank = "unknown"
		end
		write(output_file_handle, "$(rank)")
		write(output_file_handle, "\t")
		
		taxPath = map(x->split(x,"_")[3],split(unique_taxa_name,"|")); #Tax ID's
		taxPathSN = map(x->join(split(x,"_")[4:end],"_"),split(unique_taxa_name,"|")); #Taxa names
		
		#If a Tax ID is repeated at a lower taxonomic rank, this means that that rank is missing, so let's just delete it.
		for i=1:length(taxPath)
			if i>=2
				if taxPath[i] == taxPath[i-1]
					taxPath[i] = ""
					taxPathSN[i] = ""
				end
			end
		end
		#Join back up the paths
		taxPath = join(taxPath,"|")
		taxPathSN = join(taxPathSN,"|")
		
		write(output_file_handle, "$(taxPath)")
		write(output_file_handle, "\t")
		

		write(output_file_handle, "$(taxPathSN)")
		write(output_file_handle, "\t")
		write(output_file_handle, "$(taxa_abundances[unique_taxa_name])")
		write(output_file_handle, "\n")
	end
end

#Close the output file
close(output_file_handle)





