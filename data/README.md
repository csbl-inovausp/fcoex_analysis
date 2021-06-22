# Data collection

The pbmc3k dataset was directly obtained from the SeuratData package.

Gene sets (ReactomePathways.gmt) from reactome were manually downloaded from <https://reactome.org/download-data> in 2021-07-22. 

```
wget -O data/reactome.zip  https://reactome.org/download/current/ReactomePathways.gmt.zip
cd data
unzip reactome.zip
rm reactome.zip
```