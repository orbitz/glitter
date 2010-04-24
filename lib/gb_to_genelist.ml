
let gene_boundaries = Genbank.read Sys.argv.(1)

let () = Gene_predictor.write_gene_list gene_boundaries stdout
