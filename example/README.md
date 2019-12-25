# DEBKS: a tool to detect differentially expressed circular RNA

## Example


### Simulate RNA-seq
```
make simu GTF=$gtfFile  REF_DIR=$chrosomeDIR
```

### Test1: DEBKS with raw RNA-seq
```
make test1 GTF=$gtfFile REF=$genomeFa STARINDEX=$starIndex
```

### Test2: DEBKS with junction counts
```
make mapping STARINDEX=$starIndex
make test1 GTF=$gtfFile REF=$genomeFa
```
