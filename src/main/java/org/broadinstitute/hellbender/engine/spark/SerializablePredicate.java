package org.broadinstitute.hellbender.engine.spark;

import java.io.Serializable;
import java.util.function.Predicate;

@FunctionalInterface
public interface SerializablePredicate<L> extends Predicate<L>, Serializable {
}
