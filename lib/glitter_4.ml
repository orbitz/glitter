
let data = Gene_prediction_4.predict Sys.argv.(2) Sys.argv.(1)

let () = Gene_predictor.print_data stdout data


