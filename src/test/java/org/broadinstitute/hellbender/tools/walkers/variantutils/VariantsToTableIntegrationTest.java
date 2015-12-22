package org.broadinstitute.hellbender.tools.walkers.variantutils;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public final class VariantsToTableIntegrationTest extends CommandLineProgramTest {
    private String variantsToTableCmd(String moreArgs) {
        return  " --variant " + getToolTestDataDir() + "soap_gatk_annotated.noChr_lines.vcf" +
                " -T VariantsToTable" +
                " -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F TRANSITION -F DP -F SB -F set -F RankSumP -F refseq.functionalClass*" +
                " -O %s" + moreArgs;
    }

    private String variantsToTableMultiAllelicCmd(String moreArgs) {
        return  " --variant " + getToolTestDataDir() + "multiallelic.vcf" +
                " -T VariantsToTable" +
                " -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F MULTI-ALLELIC -F AC -F AF" +
                " -O %s" + moreArgs;
    }

    private String variantsToTableCmdNoSamples(String moreArgs) {
        return  " --variant " + getToolTestDataDir() + "vcfexample.noSamples.vcf" +
                " -T VariantsToTable" +
                " -O %s" + moreArgs;
    }

    @Test
    public void testInputFileFail() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant does_not_exist.vcf -o %s",
                1, UserException.CouldNotReadInputFile.class);
        spec.executeTest("testComplexVariantsToTable-FAIL", this);
    }

    @Test
    public void testOutputFileFail() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --variant " + getToolTestDataDir() + "soap_gatk_annotated.noChr_lines.vcf" +
                " -O /does_not_exists/txt.table",
                1, UserException.CouldNotCreateOutputFile.class);
        spec.executeTest("testComplexVariantsToTable-FAIL", this);
    }

    @Test
    public void testComplexVariantsToTableFail() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                variantsToTableCmd(""),
                1, UserException.class);
        spec.executeTest("testComplexVariantsToTable-FAIL", this);
    }

    @Test
    public void testNoSamples() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                variantsToTableCmdNoSamples(" -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F TRANSITION -F EVENTLENGTH"),
                Arrays.asList(getToolTestDataDir() + "expected.vcfexample.noSamples.table"));
        spec.executeTest("testNoSamples", this);
    }

    @Test
    public void testNoSamplesSoNoGenotypes() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                variantsToTableCmdNoSamples(" -GF DP"),
                1,
                UserException.class);
        spec.executeTest("testNoSamples", this);
    }

    @Test
    public void testComplexVariantsToTable() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                variantsToTableCmd(" -AMD"),
                Arrays.asList(getToolTestDataDir() + "expected.soap_gatk_annotated.noChr_lines.AMD.table"));
        spec.executeTest("testComplexVariantsToTable", this);
    }

    @Test
    public void testMultiAllelicOneRecord() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                variantsToTableMultiAllelicCmd(""),
                Arrays.asList(getToolTestDataDir() + "expected.multiallelic.table"));
        spec.executeTest("testMultiAllelicOneRecord", this);
    }

    @Test
    public void testMultiAllelicSplitRecords() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                variantsToTableMultiAllelicCmd(" -SMA"),
                Arrays.asList(getToolTestDataDir() + "expected.multiallelic.SMA.table"));
        spec.executeTest("testMultiAllelicSplitRecords", this);
    }

    @Test
    public void testGenotypeFields() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " --variant " + getToolTestDataDir() + "vcfexample2.vcf" +
                        " -T VariantsToTable" +
                        " -GF RD" +
                        " -o %s",
                Arrays.asList(getToolTestDataDir() + "expected.vcfexample2.GF_RD.table"));
        spec.executeTest("testGenotypeFields", this);
    }

    @Test
    public void testMultiallelicGenotypeFields() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " --variant " + getToolTestDataDir() + "multiallelic_gt.vcf" +
                        " -T VariantsToTable" +
                        " -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F MULTI-ALLELIC" +
                        " -GF PL -GF AD" +
                        " -SMA" +
                        " -o %s",
                Arrays.asList(getToolTestDataDir() + "expected.multiallelic_gt.table"));
        spec.executeTest("testMultiallelicGenotypeFields", this);
    }

    @Test
    public void testGenotypeFieldsWithInline() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " --variant " + getToolTestDataDir() + "vcfexample2.vcf" +
                        " -T VariantsToTable" +
                        " -GF RD -GF GT -GF GQ" +
                        " -o %s",
                Arrays.asList(getToolTestDataDir() + "expected.vcfexample2.GF_RD.GF_GT.GF_GT.table"));
        spec.executeTest("testGenotypeFieldsWithInline", this);
    }

    @Test
    public void testListFields() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " --variant " + getToolTestDataDir() + "vcfexample.withMLE.vcf" +
                        " -T VariantsToTable" +
                        " -GF PL" +
                        " -o %s",
                Arrays.asList(getToolTestDataDir() + "expected.vcfexample.withMLE.GF_PL.table"));
        spec.executeTest("testGenotypeFields", this);
    }

    @Test
    public void testMoltenOutput() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " --variant " + getToolTestDataDir() + "vcfexample2.vcf" +
                        " -T VariantsToTable" +
                        " -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER" +
                        " --moltenize" +
                        " -o %s",
                Arrays.asList(getToolTestDataDir() + "expected.vcfexample2.moltenize.table"));
        spec.executeTest("testMoltenOutput", this);
    }

    @Test
    public void testMoltenOutputWithGenotypeFields() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " --variant " + getToolTestDataDir() + "vcfexample2.vcf" +
                        " -T VariantsToTable" +
                        " -GF RD" +
                        " --moltenize" +
                        " -o %s",
                Arrays.asList(getToolTestDataDir() + "expected.vcfexample2.GF_RD.moltenize.table"));
        spec.executeTest("testMoltenOutputWithGenotypeFields", this);
    }

    @Test
    public void testMoltenOutputWithMultipleAlleles() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                        " --variant " + getToolTestDataDir() + "multiallelic.vcf" +
                        " -T VariantsToTable" +
                        " -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F MULTI-ALLELIC -F AC -F AF" +
                        " --moltenize -SMA" +
                        " -o %s",
                Arrays.asList(getToolTestDataDir() + "expected.multiallelic.moltenize.SMA.table"));
        spec.executeTest("testMoltenOutputWithMultipleAlleles", this);
    }
}
