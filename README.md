# index_bam_by_read_id
python module for sorting, indexing and retrieving reads from BAMs using read 
IDs/names

## installation

The simplest way to install the 'index_bam_by_read_id' module and the 
'bam_reads_by_id' script is to use pip:


    pip install index_bam_by_read_id

    #or without root privileges 
    pip install index_bam_by_read_id --user
    
Alternatively you may clone this repository with git and use the setup.py
script:

    git clone https://github.com/gantzgraf/index_bam_by_read_id.git
    python setup.py install 
    #or without root privileges
    python setup.py install --user
    
pysam version 0.11.1 or higher is required.

## bam_reads_by_id usage

The bundled script contains three main methods, 'sort', 'index' and 'get' in
order to sort, index and retrieve records from a BAM. It can be used as follows:

    bam_reads_by_id sort BAM [options]                                                     
    bam_reads_by_id index BAM [options]    
    bam_reads_by_id get BAM [options]      

For full documentation run:

    bam_reads_by_id -h

Note that this program can output BAM, SAM and CRAM format when using the 'sort'
and 'get' methods but only BAM and SAM formats can be indexed and be used as
input for the 'get' command.

## using the module

The index_bam_by_read_id module contains the the IndexByReadId module and two
basic Exception classes (UnsortedBamError and OutFormatError). Documentation can
be viewed by running:

    pydoc index_bam_by_read_id


## author

Written by David A. Parry at the University of Edinburgh. 


## copyright and license

MIT License

Copyright (c) 2017 David A. Parry

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


