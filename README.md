# eIF4F Molecular Dynamics Simulation

## Version Information
- **NumPy:** 1.26.4  
- **GSD:** 2.9.0  
- **HOOMD-blue:** 3.5.0  

---

## Parameter Setup
All simulation parameters are defined in the file **`parentparam.yml`**.

Each parameter can take one of the following forms:
- **Boolean:** `ON` / `OFF`
- **String**
- **List:** multiple possible values

---

## Generating folder with given parameters
To generate simulation folders for all parameter combinations, run:

```bash
python param.py parentparam.yml
```
This script iterates through all parameter combinations and creates corresponding folders.
Each generated folder includes a JSON configuration file describing the selected parameters.


---

## Run simulation
To execute a simulation for a specific folder, use:
```bash
./runscript.sh "<folder_name>"
```
