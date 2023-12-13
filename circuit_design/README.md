## Running the verification program on a circuit design

To check a counter design, save your original state in a .csv file with the same format as in this [example](https://github.com/churchlab/CAD_bio/blob/master/src/python/circuit_design/examples/test_counter.csv). Then, you can run the cript with the following arguments: 
- Path_circuit (string): Path to the counter circuit to analyse 
- Save_path (string): Path to save the output figure and the register matrix  
- Verbose (bool): To print the verification process step by step or not  

Here is an example:  
```
verification_on_counter.py test_counter.csv tests/ False
```

To check a register design, follow the same recipe with those arguments (csv file example [here](https://github.com/churchlab/CAD_bio/blob/master/src/python/circuit_design/examples/test_register.csv)):   
- Path_circuit (string): Path to the counter circuit to analyse 
- Path_program (string): Path to the counter regulatory program  
- Save_path (string): Path to save the output figure and the register matrix  
- Name (string): Name of the output file
- Verbose (bool): To print the verification process step by step or not

```
run_verification_on_register.py test_register.csv ../../../notebooks_GRNs/designs/A030-thesis-count8/changes.txt tests/ register_output False

```
