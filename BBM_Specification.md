# BBM (Binary BisMap) Format Specification:

BBM is a binary format for storing a single data track for a genome. Its main restriction is that it can only store integer values in the range 0-100 (inclusive). This results from the fact it is used to store Bismap precent mappabilities, which fall in that range. However, any data track that can be converted into 101 discrete integer values could be stored in this format. The format achieves a small file size through run-length compression of the data, with a maximum run-length (encoded as one unit) of 65535 bases. Longer runs can be encoded as sequences of shorter runs. In a test on a Bismap human genome mappability dataset, this format achieves a compression ratio of approximately 7.7:1. The format is as follows.

## Format Overview:

(`[bracketed]` items will be explained below, `<angle bracketed>` items are simply markers of the beginning and end of the file, and items in `{curly braces}` contain other items and are explained in their own sections)

```<beginning of file>[version][chromCount]{chrom record}{chrom record}{chrom record}{chrom record}...<end of file>```

* version: an unsigned character denoting the version of the BBM format used by this file. The version described here has an ID of 1.

* chromCount: an unsigned little-endian 32-bit integer, reprensenting the number of chromosomes in the file

* chrom record: This stores the data for a given chromosome and will be discussed below.

## Chromosome Records:

A chromosome record stores the name of a chromosome and the values in the data track for that chromosome. It is structured as follows:

[nameLength][chromName][NUL][chromLength]{chrom data}

* nameLength: an unsigned little-endian 16-bit integer, reprensenting the number of characters in the chromosome's name (not including the null terminator)

* chromName: a sequence of 8-bit ASCII characters representing the name of the chromosome

* NUL: a \0 byte functioning as a null terminator for the name. The presence of this null terminator can be used for error checking.

* chromLength: an unsigned little-endian 32-bit integer, reprensenting the number of bases in the chromosome

* chrom data: This stores the values in the data track for the chromosome and will be discussed below.

## Chromosome Data:

The BBM format stores run-length-compressed data for a single data track for each chromosome. Each value (not including runs) is represented by a value of 0-100 (inclusive) in a single byte. Runs are represented in two formats: short runs (lengths from 2-154 inclusive), which are stored in two bytes; and long runs (lengths from 155-65535 inclusive), which are stored in four bytes. The value 255 is a marker for a long run and is not to be used as an individual value or a length marker for a short run.

The format of this is simply a series of individual values, short runs, and long runs. Positions are not stored explicitly, the first value is at position 0 in the chromosome and successive values are at consecutive positions, ending with the last value at position chromLength-1. There is no specific end marker for the chromosome, it ends after chromLength values have been listed.

### Individual Value

This is stored as a single unsigned byte containing a value in the range of 0-100 (inclusive), representing the value of the data track at that position.

### Short Run (length 2-155)

This is represented as follows: 

```[runLength][runValue]```

* runLength: A single unsigned byte storing a value in the range of 102-154 (inclusive), representing the run length plus 99, (i.e. 101 represents a run of length 2, 102 a run of 3, and so on until 254 (which repesents a run of 155)). As stated above, this must not contain the value 255, for reasons that will be detailed in the section on long runs.
* runValue: A single unsigned byte storing the run's value.

### Long Run (length 155-65535)

This is represented as four bytes as follows:

```[255][runLength][runValue]```

 * 255: This is a flag indicating this is a long run. This is the reason that 255 cannot be used in short runs, since its presence here at the start of a run's data indicates a long run.
 * runLength: This is a little-endian 16-bit unsigned integer representing the length of the run. 
 * runValue: An unsigned byte containing the run's value.
