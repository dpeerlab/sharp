version 1.0

import "modules/Count.wdl" as module

workflow Count {

    input {
        File fastqR1
        File fastqR2
        File cbWhiteList
        File hashTagList

        # cellular barcode
        Int cbStartPos
        Int cbEndPos

        # UMI
        Int umiStartPos
        Int umiEndPos

        Int cbCollapsingDistance
        Int umiCollapsingDistance

        Int numExpectedCells
    }

    call module.CiteSeqCount {
        input:
            fastqR1 = fastqR1,
            fastqR2 = fastqR2,
            cbWhiteList = cbWhiteList,
            hashTagList = hashTagList,
            cbStartPos = cbStartPos,
            cbEndPos = cbEndPos,
            umiStartPos = umiStartPos,
            umiEndPos = umiEndPos,
            cbCollapsingDistance = cbCollapsingDistance,
            umiCollapsingDistance = umiCollapsingDistance,
            numExpectedCells = numExpectedCells
    }

    output {
        File outUmiDenseCount = CiteSeqCount.outUmiDenseCount
        File outUnmapped = CiteSeqCount.outUnmapped
        File outReport = CiteSeqCount.outReport
        File outUncorrected = CiteSeqCount.outUncorrected

        Array[File] outUmiCount = CiteSeqCount.outUmiCount
        Array[File] outReadCount = CiteSeqCount.outReadCount
    }
}
