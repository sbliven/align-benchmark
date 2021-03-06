package org.rcsb.alignBenchmark.metrics;


import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.model.AFPChain;
import org.rcsb.alignBenchmark.metrics.Metric;
import org.rcsb.alignBenchmark.MultipleAlignment;

public abstract class AlignmentLengthMetric extends Metric{
	/**
	 * Calculates the length of the reference alignment
	 * @author Spencer Bliven
	 *
	 */
	public static class Reference extends AlignmentLengthMetric {

		@Override
		public double calculate(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {
			return reference.size();
		}

		@Override
		public String getName() {
			return "Ref_len";
		}


	}
	
	/**
	 * Calculates the length of the test alignment
	 * @author Spencer Bliven
	 *
	 */
	public static class Alignment extends AlignmentLengthMetric {

		@Override
		public double calculate(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {
			return (double)align.getOptLength();
		}

		@Override
		public String getName() {
			return "Aln_len";
		}
	}
	
	@Override
	public String format(double result) {
		return Integer.toString((int)result);
	}
}
