# ./run.py --method [GSEA|iPanda|iF|conexic] --control control.ext --case tumor.ext
# ./run.py --method [GSEA|iPanda|iF|conexic] --control some_file.ext --case 20degree.ext 30degree.ext
# ./run.py --method [GSEA|iPanda|iF|conexic] --data some_file.ext --case :4    # reszta to kontrole
# ./run.py --method [GSEA|iPanda|iF|conexic] --data some_file.ext --control :3
# ./run.py --method [GSEA|iPanda|iF|conexic] --data some_file.ext --control 1,3 # all other are cases
# ./run.py --method [GSEA|iPanda|iF|conexic] --data some_file.ext --case 1,2,4 --control 5-8 # olewamy trójkę
# ./run.py --method [GSEA|iPanda|iF|conexic] --data some_file.ext --case 1,2,4 --control 5:8

# ./run.py --method [GSEA|iPanda|iF|conexic] --control controls.ext --columns 1:4 --case case.ext --columns 1:3

# TODO: czy większość danych ma headery? Czy switch powinien być domyślnie włączony?
# ./run.py --method [GSEA|iPanda|iF|conexic] --control controls.ext --name '20stopni' --columns 1:4 --do_not_use_sample_names_from_header\
#    --case case.ext --name '40stopni' --columns 1:3
# ./run.py --method [GSEA|iPanda|iF|conexic] --control controls.ext 1:4 --case case.ext 1-3
# ./run.py --method [GSEA|iPanda|iF|conexic] --control controls.ext 1:4 --case '20stopni' case.ext 1-3 --case '40stopni' case.ext 4-5

# ./run.py --method iPanda --some_method_argument some_value --control
