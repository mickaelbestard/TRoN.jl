### itérations 0 à 15: 
    -   descente de gradient lente, on sort du minimum local grâce à FP.

### itérations 15 à 30:
    -   GD rapide, on est peut-être dans un bon bassin d'attraction.

### itérations 30 à 42(fin):
    -   FP avait gêné la convergence, mais on est restés dans un bassin
        qui permet une convergence rapide vers une loss basse.



### Remarque: il peut être judicieux de déclencher FP de façon astucieuse :
    -   toutes les K itérations, puis 2K, 3K...
    -   en fonction de la valeur de la loss (se fixer un seuil ɛ tel que le bassin
        d'attraction est considéré comme "bon" pour loss < ɛ, et ne plus utiliser 
        le point fixe à partir de là).