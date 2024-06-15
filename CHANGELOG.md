# Changelog

## [3.0.3](https://github.com/cuspuk/workflow_mapping/compare/v3.0.2...v3.0.3) (2024-06-15)


### Performance Improvements

* autobump conda envs and/or wrappers ([171312f](https://github.com/cuspuk/workflow_mapping/commit/171312fc760a9f7c737a84629e41bf82dad2c520))

## [3.0.2](https://github.com/cuspuk/workflow_mapping/compare/v3.0.1...v3.0.2) (2024-06-08)


### Performance Improvements

* autobump conda envs and/or wrappers ([e6564ae](https://github.com/cuspuk/workflow_mapping/commit/e6564ae0d234d64b3559b0f485b8e763018716d6))

## [3.0.1](https://github.com/cuspuk/workflow_mapping/compare/v3.0.0...v3.0.1) (2024-05-15)


### Bug Fixes

* changes after repo transfer ([e672167](https://github.com/cuspuk/workflow_mapping/commit/e67216766139653904c7d00bb031c967a7e2197a))

## [3.0.0](https://github.com/xsitarcik/mapping/compare/v2.3.2...v3.0.0) (2024-05-10)


### ⚠ BREAKING CHANGES

* multiqc is for each reference separately

### Features

* multiqc is for each reference separately ([1993edf](https://github.com/xsitarcik/mapping/commit/1993edfb278b5d626385a42ddd6760244917ee10))


### Bug Fixes

* added multiqc into report ([0f0e7b1](https://github.com/xsitarcik/mapping/commit/0f0e7b10365d68f1225207c4d4eb9d5917039da4))
* multiqc updates ([65f7c7b](https://github.com/xsitarcik/mapping/commit/65f7c7b6a6924ce68cb3e1dabb39bfa29b4109c6))
* reads module output added into get_outputs ([7f4ddf9](https://github.com/xsitarcik/mapping/commit/7f4ddf96bef8e99136ce204a5024cff61acfefff))
* restructured outputs from step/sample to sample.step ([0df7f63](https://github.com/xsitarcik/mapping/commit/0df7f6392164e95d4542610ac9ea609526949545))
* restructured outputs from step/sample to sample.step ([0834f60](https://github.com/xsitarcik/mapping/commit/0834f601e402d0de3a77e4e0cd805f958329a0ea))


### Performance Improvements

* bumped wrappers ([682a375](https://github.com/xsitarcik/mapping/commit/682a375bbc5d3998aa209a2508492e5b6adff113))

## [2.3.2](https://github.com/xsitarcik/mapping/compare/v2.3.1...v2.3.2) (2024-04-30)


### Bug Fixes

* multiqc rule takes config as input ([167a228](https://github.com/xsitarcik/mapping/commit/167a22812043cd8c47afe91b6f0f0b152995bb68))
* report reorganization ([1bef590](https://github.com/xsitarcik/mapping/commit/1bef590d3d98080c73c3e9b91ef5f87e070562ac))
* small changes in schema to prevent bad inputs for resources ([824bb55](https://github.com/xsitarcik/mapping/commit/824bb55fa7b46cdd9f8edb1af8a24584f40a8722))


### Performance Improvements

* bumped reads module ([2aadfe8](https://github.com/xsitarcik/mapping/commit/2aadfe8b9622ff511504d5904c6bdab6cb0d1ac6))
* bumped wrappers ([26f2cd0](https://github.com/xsitarcik/mapping/commit/26f2cd0fffc09621e8a1a09a8b19b1cdff222030))

## [2.3.1](https://github.com/xsitarcik/mapping/compare/v2.3.0...v2.3.1) (2024-04-25)


### Bug Fixes

* mapped bam is not temporary if deduplication is not used ([6b4c149](https://github.com/xsitarcik/mapping/commit/6b4c1498f01a727b4707f04c6ffbf6ff0ea2b3f8))

## [2.3.0](https://github.com/xsitarcik/mapping/compare/v2.2.0...v2.3.0) (2024-04-25)


### Features

* added samtools stats ([5d55bbe](https://github.com/xsitarcik/mapping/commit/5d55bbe05797e289a6f21611db40935595fe26e2))
* multiqc ([799e3c4](https://github.com/xsitarcik/mapping/commit/799e3c499558f67d4acfd6f15e8b37c50b7c689d))


### Bug Fixes

* qualimap produced into bamqc/sample instead of reverse ([f6e39e6](https://github.com/xsitarcik/mapping/commit/f6e39e6ed020195aa7679026816ff14433611114))

## [2.2.0](https://github.com/xsitarcik/mapping/compare/v2.1.0...v2.2.0) (2024-04-22)


### Features

* fasta paths required instead of dirs for references ([a0e3e99](https://github.com/xsitarcik/mapping/commit/a0e3e99e288b58a187e6854d7f7b28ae7a5c196d))


### Performance Improvements

* bumped wrappers ([8609298](https://github.com/xsitarcik/mapping/commit/8609298b46a356c3662049ea911040842ad79a4c))

## [2.1.0](https://github.com/xsitarcik/mapping/compare/v2.0.2...v2.1.0) (2024-04-16)


### Features

* added module versions configurable ([7ebf26e](https://github.com/xsitarcik/mapping/commit/7ebf26e01d0313e37e8bb870428fb17aa74873a1))

## [2.0.2](https://github.com/xsitarcik/mapping/compare/v2.0.1...v2.0.2) (2024-02-25)


### Bug Fixes

* updated schemas ([180c5eb](https://github.com/xsitarcik/mapping/commit/180c5eb595ad981e96e01514ad34eba216d4dbf0))

## [2.0.1](https://github.com/xsitarcik/mapping/compare/v2.0.0...v2.0.1) (2024-01-26)


### Bug Fixes

* qualimap requested for steps defined in config ([99c558a](https://github.com/xsitarcik/mapping/commit/99c558af84978c9307a47a8d59f8694388b393a0))

## [2.0.0](https://github.com/xsitarcik/mapping/compare/v1.0.0...v2.0.0) (2024-01-26)


### ⚠ BREAKING CHANGES

* added method specifics in config

### Features

* added method specifics in config ([c7c1caa](https://github.com/xsitarcik/mapping/commit/c7c1caadd3444997e7fd2bad2284e5648ae5b1cf))
* expanded schema to dynamically restrict config ([3945b13](https://github.com/xsitarcik/mapping/commit/3945b137e648ec053b192a017b6a40009a30bf0e))

## 1.0.0 (2024-01-05)


### Features

* first version of mapping ([d9520ab](https://github.com/xsitarcik/mapping/commit/d9520ab605bc538a4a9b1b3d4236e26641c64f28))
