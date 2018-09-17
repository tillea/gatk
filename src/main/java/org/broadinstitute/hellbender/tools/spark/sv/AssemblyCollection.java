package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.TemplateFragmentOrdinal;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.WeakHashMap;

/**
 * Represents a collection of assembly files  of file in the assembly collection
 */
public class AssemblyCollection implements Serializable {

    private static final long serialVersionUID = 1L;

    private final String location;
    private final String fastqFileFormat;

    private final WeakHashMap<Integer, List<Template>> templatesById;

    public AssemblyCollection(final String location, final String fastqFileFormat) {
        this.location = Utils.nonNull(location);
        this.fastqFileFormat = Utils.nonNull(fastqFileFormat);
        this.templatesById = new WeakHashMap<>(10);
    }

    public List<Template> templates(final int assemblyNumber) {
        ParamUtils.isPositiveOrZero(assemblyNumber, "the assembly number must be 0 or greater");
        // Is important to create a separate key Integer object so that WeakHash entry can be re-claimed.
        // Auto-boxing may reuse same instance for numbers close enough to 0 (see {@link Integer#valueOf} for details).
        @SuppressWarnings("UnnecessaryBoxing") final Integer key = new Integer(assemblyNumber);
        return templatesById.computeIfAbsent(assemblyNumber, this::readTemplates);
    }

    private List<Template> readTemplates(final int assemblyNumber) {
        Utils.nonNull(assemblyNumber);
        final String path = IOUtils.getPath(location).resolve(String.format(fastqFileFormat, assemblyNumber)).toString();
        final List<Template> result = new ArrayList<>();
        if (!BucketUtils.fileExists(path)) {
            throw new UserException.CouldNotReadInputFile(path, "missing input file for assembly number " + assemblyNumber);
        } else {
            final List<SVFastqUtils.FastqRead> fastqReads = SVFastqUtils.readFastqFile(path.toString());
            final int numberOfFastqReads = fastqReads.size();
            if ((numberOfFastqReads & 1) == 1) {
                throw new UserException.BadInput("There is a odd number of reads in the assembly fastq file: " + path);
            }
            for (int i = 0; i < fastqReads.size(); i += 2) {
                final SVFastqUtils.FastqRead first = fastqReads.get(i);
                final SVFastqUtils.FastqRead second = fastqReads.get(i + 1);
                if (!first.getName().equals(second.getName())) {
                    throw new UserException.BadInput(
                            String.format("consecutive reads don't have the same name '%s' != '%s' in '%s'", first.getName(), second.getName(), path));
                } else if (first.getFragmentOrdinal() != TemplateFragmentOrdinal.PAIRED_FIRST) {
                    throw new UserException.BadInput(
                            String.format("first in pair fragment number is not 1 (%s) for '%s' in '%s'", first.getFragmentOrdinal().nameSuffix(), first.getName(), path));
                } else if (second.getFragmentOrdinal() != TemplateFragmentOrdinal.PAIRED_SECOND) {
                    throw new UserException.BadInput(
                            String.format("first in pair fragment number is not 2 (%s) for '%s' in '%s'", second.getFragmentOrdinal().nameSuffix(), first.getName(), path));
                } else {
                    final Template template = Template.create(first.getName(), Arrays.asList(first, second),
                            Template.Fragment::of);
                    result.add(template);
                }
            }
        }
        return result;
    }
}
