
let ((total, path, prob), gene_list, count, training_count) = Gene_prediction_10.predict Sys.argv.(2) Sys.argv.(1)

let () =
  Printf.fprintf stdout "%f\n%f\n%d\n%d\n" total prob count training_count;
  Gene_predictor.write_gene_list gene_list stdout


