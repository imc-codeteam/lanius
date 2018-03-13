# Lanius

Lanius is a small program to build crystal slab models for VASP.

## Usage

To construct a six-layer fcc(111) surface of Co atoms with a vacuum size of 10 Angstrom:

```
./lanius -m 1,1,1 -d 3,3,6 -e Co -c fcc -o POSCAR
```

### Input parameters

* `-m` Miller index, e.g. `1,1,1`
* `-d` Unit cell dimensions, e.g. `2,2,4`
* `-e` Element, e.g. `Co`
* `-v` Vacuum thickness in Angstrom, e.g. `10`
* `-c` Crystal structure, choose from `bcc`, `hcp` or `fcc`
* `-m` Miller index, e.g. `1,1,1`

## Compilation

```
mkdir build
cd build
cmake ../src
make -j4
```
