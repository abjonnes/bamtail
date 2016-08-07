bamtail
=======

A simple tool to quickly tell you about the final-ish alignment present in a [BAM file](http://genome.sph.umich.edu/wiki/BAM) at runtime. This is useful for monitoring progress as a large, position-sorted BAM file is being generated. Since BAM files are compressed, ``tail`` is not useful for this! ``bamtail`` is compatible with both ``python2`` and ``python3``.


Usage
-----

```
$ bamtail sample1.bam
chr4:34062067

$ bamtail sample2.bam sample3.bam sample4.bam
sample2.bam: chr3:23469490
sample3.bam: unmapped
sample4.bam: chr14:1045735
```


Installation
------------

You can install with `pip`:
```
pip install https://github.com/abjonnes/bamtail/tarball/master
```
Or, with `setuptools` installed:
```
git clone https://github.com/abjonnes/bamtail.git
cd bamtail
python setup.py install
```


Changelog
---------
- compatibility with py2
