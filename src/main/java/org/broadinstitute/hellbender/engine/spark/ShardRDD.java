package org.broadinstitute.hellbender.engine.spark;

import com.google.common.collect.Iterators;
import org.apache.spark.HashPartitioner;
import org.apache.spark.Partitioner;
import org.apache.spark.api.java.AbstractJavaRDDLike;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaRDDLike;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.engine.ShardBoundaryShard;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SerializableSupplier;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;
import shapeless.ops.hlist;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * {@link JavaPairRDD} of Shards where the key is the shard interval and the value is
 * the record of interest at each shard.
 * @param <L>
 */
public class ShardRDD<L> implements Serializable {

    private static final long serialVersionUID = 1L;

    private JavaPairRDD<ShardBoundary, List<L>> pairRDD;

    public JavaPairRDD<ShardBoundary, List<L>> toPairRDD() {
        return pairRDD;
    }

    /**
     * Removed duplicated record within each {@link Shard} record
     * list.
     * <p>
     *     It relies on record's implementation of {@link Object#hashCode}
     *     and {@link Object#equals}.
     * </p>
     *
     * @return never {@code null}.
     */
    public ShardRDD<L> distinctWithinShard() {
        return new ShardRDD<>(pairRDD.mapPartitionsToPair(it -> {
            final LinkedHashSet<L> unique = new LinkedHashSet<>();
            return Utils.map(it, tuple -> {
                unique.clear();
                boolean repeatsFound = false;
                for (final L record : tuple._2) {
                    if (!unique.add(record)) {
                        repeatsFound = true;
                    }
                }
                // if no repeats are found we can simply use the original
                // tuple.
                return !repeatsFound ? tuple :
                        new Tuple2<>(tuple._1, new ArrayList<>(unique));
            });
        }, true));
    }

    /**
     * Returns a new ShardRDD that merges to record of this and another
     * RDD within each shard.
     * <p>
     *     Notice that for the join to work properly both RDD should share
     *     the same shards. For shard records to be merged the shard boundaries have to
     *     be exactly the same in both RDDs.
     * </p>
     *
     * @param other
     * @return never {@code null}.
     */
    public ShardRDD<L> join(final ShardRDD<L> other) {
        Utils.nonNull(other);
        return new ShardRDD<>(
                pairRDD.cogroup(other.pairRDD)
                .mapValues(tuple ->
                        Stream.concat(Utils.stream(tuple._1)
                                .flatMap(List::stream),
                                Utils.stream(tuple._2).flatMap(List::stream))
                                .collect(Collectors.toList())));
    }

    public interface Grouper<V, W> extends Serializable {
        List<Tuple2<V, List<W>>> group(final List<V> left, final List<W> right);
    }

    ShardRDD(final JavaPairRDD<ShardBoundary, List<L>> pairRDD) {
        this.pairRDD = Utils.nonNull(pairRDD);
    }

    public <W> ShardRDD<W> mapRecords(final SerializableFunction<L, W> mapper) {
        return new ShardRDD<>(pairRDD.mapValues(list -> list.stream().map(mapper).collect(Collectors.toList())));
    }

    public <W> ShardRDD<W> flatMapRecords(final SerializableFunction<L, Stream<W>> mapper) {
       return new ShardRDD<>(pairRDD.mapValues(list ->
           list.stream().flatMap(mapper).collect(Collectors.toList())));
    }

    public ShardRDD<L> filterRecords(final SerializablePredicate<L> predicate) {
        return new ShardRDD<>(pairRDD.mapPartitionsToPair(inShardTuplesIt ->
            Utils.map(inShardTuplesIt, tuple -> {
                final List<L> oldRecords = tuple._2;
                final List<L> newRecords = oldRecords.stream().filter(predicate).collect(Collectors.toList());
                if (oldRecords.size() == newRecords.size()) {
                    return tuple;
                } else {
                    return new Tuple2<>(tuple._1, newRecords);
                }
            }), true));
    }

    /**
     * Transforms the shard records by mapping the entire record list per shard while keeping the shard boundaries unchanged.
     * @param mapperSupplier methods that will supplier the worker mapper function.
     * @param <W> the type-parameter of the output shard record.
     * @return never {@code null}.
     */
    public <W> ShardRDD<W> mapPartitionShardRecords(final SerializableSupplier<BiFunction<ShardBoundary, List<L>, List<W>>> mapperSupplier) {
        return new ShardRDD<>(pairRDD.mapPartitionsToPair(inShardTuplesIt -> {
            final BiFunction<ShardBoundary, List<L>, List<W>> mapper = mapperSupplier.get();
            return Utils.map(inShardTuplesIt, tuple -> new Tuple2<>(tuple._1, mapper.apply(tuple._1, tuple._2)));
        }, true));
    }

    /**
     * Transforms the shard records mapping them one by one while keeping the shard boundaries unchanged.
     * <p>
     *     The mapper supplied is a second order lambda function that composes and returns the actual mapper function.
     *     This is needed in order to enable the caller to incorporate once-per-partition initialization code
     *     while preventing the caller code to change the shard boundaries (keys).
     * </p>
     * @param mapperSupplier methods that will supplier the worker mapper function.
     * @param <W> the type-parameter of the output shard record.
     * @return never {@code null}.
     */
    public <W> ShardRDD<W> mapPartitionRecords(final SerializableSupplier<Function<L, W>> mapperSupplier) {
        return new ShardRDD<>(pairRDD.mapPartitionsToPair(inShardTuplesIt -> {
            final Function<L, W> mapper = mapperSupplier.get();
            return Utils.map(inShardTuplesIt, tuple -> new Tuple2<>(tuple._1, tuple._2().stream().map(mapper).collect(Collectors.toList())));
        }, true));
    }

    public <W> ShardRDD<Tuple2<L, List<W>>> groupRight(final ShardRDD<W> right, final Grouper<L, W> grouper) {
        return new ShardRDD<>(pairRDD.join(right.pairRDD).mapValues(values -> grouper.group(values._1(), values._2())));
    }

    /**
     * Returns an {@link JavaRDD} of all the records flatten into a single stream.
     * @return never {@code null}.
     */
    public JavaRDD<L> records() {
        return pairRDD.values().flatMap(List::iterator);
    }
}
