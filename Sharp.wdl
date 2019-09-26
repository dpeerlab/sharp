version 1.0

import "modules/MergeFastq.wdl" as MergeFastq
import "modules/FastQC.wdl" as FastQC
import "modules/Cutadapt.wdl" as Cutadapt
import "modules/PrepCBWhitelist.wdl" as PrepCBWhitelist
import "modules/Count.wdl" as Count
import "modules/HtoDemux.wdl" as HtoDemux
import "modules/Combine.wdl" as Combine

workflow Sharp {

    input {
        Array[File] uriFastqR1
        Array[File] uriFastqR2

        Int lengthR1
        Int lengthR2

        String sampleName

        File cellBarcodeWhitelistUri
        String cellBarcodeWhiteListMethod

        File hashTagList

        # cellular barcode start/end positions
        Int cbStartPos
        Int cbEndPos

        # UMI start/end positions
        Int umiStartPos
        Int umiEndPos

        # correction
        Int cbCollapsingDistance
        Int umiCollapsingDistance

        Int numExpectedCells

        Int numCoresForCount

        Float htoDemuxQuantile

        File denseCountMatrix
    }

    # merge FASTQ R1
    call MergeFastq.MergeFastq as MergeFastqR1 {
        input:
            uriFastq = uriFastqR1,
            sampleName = sampleName,
            readName = "R1"
    }

    # merge FASTQ R2
    call MergeFastq.MergeFastq as MergeFastqR2 {
        input:
            uriFastq = uriFastqR2,
            sampleName = sampleName,
            readName = "R2"
    }

    # run FastQC R1
    call FastQC.FastQC as FastQCR1 {
        input:
            fastqFile = MergeFastqR1.out
    }

    # run FastQC R2
    call FastQC.FastQC as FastQCR2 {
        input:
            fastqFile = MergeFastqR2.out
    }

    # trim R1
    call Cutadapt.Trim as TrimR1 {
        input:
            fastq = MergeFastqR1.out,
            length = lengthR1,
            outFileName = "R1.fastq.gz"
    }

    # trim R2
    call Cutadapt.Trim as TrimR2 {
        input:
            fastq = MergeFastqR2.out,
            length = lengthR2,
            outFileName = "R2.fastq.gz"
    }

    # prepare cell barcode whitelist
    if (cellBarcodeWhiteListMethod == "SeqcSparseCountsBarcodesCsv") {
        # *_sparse_counts_barcodes.csv
        call PrepCBWhitelist.WhitelistFromSeqcSparseBarcodes {
            input:
                csvFile = cellBarcodeWhitelistUri
        }
    }

    if (cellBarcodeWhiteListMethod == "SeqcDenseCountsMatrixCsv") {
        # *_dense.csv
        call PrepCBWhitelist.WhitelistFromSeqcDenseMatrix {
            input:
                csvFile = cellBarcodeWhitelistUri
        }
    }

    if (cellBarcodeWhiteListMethod == "BarcodeWhitelistCsv") {
        # one barcode per line
        call PrepCBWhitelist.NotImplemented
    }

    File cbWhitelist = select_first([WhitelistFromSeqcSparseBarcodes.out, WhitelistFromSeqcDenseMatrix.out])

    # run CITE-seq-Count
    call Count.CiteSeqCount {
        input:
            fastqR1 = TrimR1.outFile,
            fastqR2 = TrimR2.outFile,
            cbWhiteList = cbWhitelist,
            hashTagList = hashTagList,
            cbStartPos = cbStartPos,
            cbEndPos = cbEndPos,
            umiStartPos = umiStartPos,
            umiEndPos = umiEndPos,
            cbCollapsingDistance = cbCollapsingDistance,
            umiCollapsingDistance = umiCollapsingDistance,
            numExpectedCells = numExpectedCells,
            numCores = numCoresForCount
    }

    # HTO demux
    call HtoDemux.HtoDemux {
        input:
            umiCountFiles = CiteSeqCount.outUmiCount,
            quantile = htoDemuxQuantile
    }

    # combine count matrix with hashtag
    call Combine.HashedCountMatrix {
        input:
            denseCountMatrix = denseCountMatrix,
            htoDemuxMatrix = HtoDemux.outClassCsv
    }

    # correct false positive doublets
    call Combine.CorrectFalsePositiveDoublets {
        input:
            htoClassification = HashedCountMatrix.outClass,
            denseCountMatrix = denseCountMatrix,
            umiCountFiles = CiteSeqCount.outUmiCount
    }


    output {
        File fastQCR1Html = FastQCR1.outHtml
        File fastQCR2Html = FastQCR2.outHtml

        File countReport = CiteSeqCount.outReport

        File htoClassification = HashedCountMatrix.outClass
        File hashedCountMatrix = HashedCountMatrix.outCountMatrix
        File statsClassification = HashedCountMatrix.outStats
        File logCombine = HashedCountMatrix.outLog

        File htoCorrectedClassification = CorrectFalsePositiveDoublets.outClass
        File hashedCorrectedCountMatrix = CorrectFalsePositiveDoublets.outCountMatrix
        File statsCorrectedClassification = CorrectFalsePositiveDoublets.outStats
        File logFPDoubletCorrection = CorrectFalsePositiveDoublets.outLog
    }
}
