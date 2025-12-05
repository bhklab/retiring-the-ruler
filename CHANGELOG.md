# Changelog

## [0.2.0](https://github.com/bhklab/retiring-the-ruler/compare/v0.1.0...v0.2.0) (2025-12-05)


### Features

* add assessment for different target selections, saving out synthetic data ([6331710](https://github.com/bhklab/retiring-the-ruler/commit/633171061222eeb4006c3c4f35e36816f485b62d))
* add check for non-empty dataframes in recist metrics function ([0f30fd5](https://github.com/bhklab/retiring-the-ruler/commit/0f30fd591ea5d38b982c8188291e899fd5290164))
* add click inputs ([c910f9c](https://github.com/bhklab/retiring-the-ruler/commit/c910f9cc8687fb91d6a2e0012cf6aab4e5954106))
* add logging to main pipeline code ([7958340](https://github.com/bhklab/retiring-the-ruler/commit/7958340d2f6def03a740681ee567dcf0d1d167cb))
* add parallel lesion generation ([29a3473](https://github.com/bhklab/retiring-the-ruler/commit/29a34731350db4362386c0e73b300a2e48b1b772))
* add parallel processing to select target lesions ([590be1f](https://github.com/bhklab/retiring-the-ruler/commit/590be1f9bd48c72e26bafb6ad68dc01556d3aa3e))
* add pd sensitivity calculation to recist metrics function ([7b3a394](https://github.com/bhklab/retiring-the-ruler/commit/7b3a3940ff0061842e1e22840b0a42d67836b610))
* add plotting functions ([9a394ae](https://github.com/bhklab/retiring-the-ruler/commit/9a394aeac090bb4347a39d71b39db07a84820998))
* add recist assessment functionality ([f3d91a2](https://github.com/bhklab/retiring-the-ruler/commit/f3d91a2435f71619e348db446433da17a18e7f5c))
* add volume calculations to synthetic lesion generation ([3355304](https://github.com/bhklab/retiring-the-ruler/commit/33553045e67f826690911dcde885a0d6d20ce51e))
* add volume vs diameter plots ([9f55068](https://github.com/bhklab/retiring-the-ruler/commit/9f5506871fe6905ec396ba92bd33a80a01726e91))
* function for selecting target lesions and calculating accuracy of response by lesion count ([a05cc70](https://github.com/bhklab/retiring-the-ruler/commit/a05cc7063fda1fdeec21a9bd71b59a27b766d937))
* make default value for random seed None in click ([791f4fb](https://github.com/bhklab/retiring-the-ruler/commit/791f4fb144e3dd2c39b15b10450b6b7a4848077a))
* make pixi task for seed and no seed pipeline runs ([67558f6](https://github.com/bhklab/retiring-the-ruler/commit/67558f6dbd3fa0d8295197a509c74cedc9b508e8))
* parallelize recist assessment ([0d5d704](https://github.com/bhklab/retiring-the-ruler/commit/0d5d7046239e50b03a054da8cbdf03bb5dcdd32e))
* reconfigure diameter change generation to be reproducible ([ee029d0](https://github.com/bhklab/retiring-the-ruler/commit/ee029d02f1cdbabed9a16af65a8f58f9e72dbeaa))
* set random seeds for all random choices, update random to rng ([47ef40e](https://github.com/bhklab/retiring-the-ruler/commit/47ef40e15178f6c6a7dceeacd5a31777d326cac7))
* start pipeline code, generate synthetic lesions ([91b3a41](https://github.com/bhklab/retiring-the-ruler/commit/91b3a4176a6a676424045e904f8986517b513619))


### Bug Fixes

* correct n_lesion selection so value is between 1 and 30, add patient id and lesion index to lesion data return ([b5f1d3f](https://github.com/bhklab/retiring-the-ruler/commit/b5f1d3fe2413dc9b37aaf139f956bed8b10b39d4))
* correct where the random number generator should have a random seed and pass the same rng for lesion selection tasks, intialize only once ([bf3df8b](https://github.com/bhklab/retiring-the-ruler/commit/bf3df8bed6d5bcdc1be65433dee4178aa88efa70))
* leave patient_id as number for sorting purposes ([5f23b79](https://github.com/bhklab/retiring-the-ruler/commit/5f23b7908839296712e01149e8d910ee52508718))
* use rng initialized outside of function for poisson selection ([9daaa91](https://github.com/bhklab/retiring-the-ruler/commit/9daaa91fe5afef5957ebb4028fc55a09ae3b43b8))
