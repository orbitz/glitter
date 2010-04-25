
let data = Gene_prediction_7.predict Sys.argv.(2) Sys.argv.(1)

let () = Gene_predictor.print_data stdout data


