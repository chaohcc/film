# film
a fully learned index for larger-than-memory databases


## the preformance of record size
the model size comparison of differnet methods in terms of record size.

<img src="https://user-images.githubusercontent.com/51820918/155710404-2abb0e9a-7a74-4718-b47c-8bceddbb463c.png" width="50%" height="50%">


![recordSizewiki_ts_add](https://user-images.githubusercontent.com/51820918/155705150-5a7aa409-503d-4ef0-9e06-ef00f2fc7db8.png)

## the range query performance with different amount of available memory 

<img src="https://user-images.githubusercontent.com/51820918/155710245-68bd16c0-8e0c-487d-9b74-51d34bd0871b.png" width="50%" height="50%">


## dataset
the books and wiki_ts are come from SOSD. ref: https://github.com/learnedsystems/SOSD

the optimal solution of generating piece-wise-linear functions has well studied by computional geometry [ref: Joseph O’Rourke. 1981. An on-line algorithm for fitting straight lines between data ranges. Commun. ACM 24, 9 (1981), 574–578.], and the PGM-Index has implemented it in C++ implementation[ref: https://github.com/gvinciguerra/PGM-index], the learned model of FILM is on the basis of the segmentation from PGM-index.

