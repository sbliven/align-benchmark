package org.rcsb.alignBenchmark.metrics;


import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.SVDSuperimposer;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.jama.Matrix;
import org.rcsb.alignBenchmark.MultipleAlignment;


public abstract class TMScoreMetric extends Metric{
   public static class Reference extends TMScoreMetric {

      @Override
      public double calculate(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {
         try {

            List<Atom[]> structures = new ArrayList<Atom[]>(2);
            structures.add(ca1);
            structures.add(ca2);
            int[][] optAln = reference.getAlignmentMatrix(structures);

            // Create new arrays for the subset of atoms in the alignment.
            Atom[] ca1aligned = new Atom[reference.size()];
            Atom[] ca2aligned = new Atom[reference.size()];
            int pos=0;
            for(int i=0;i<optAln[0].length;i++) {
               ca1aligned[pos] = ca1[optAln[0][pos]];
               ca2aligned[pos] = (Atom) ca2[optAln[1][pos]].clone();
               pos++;

            }

            //Superimpose
            SVDSuperimposer svd = new SVDSuperimposer(ca1aligned, ca2aligned);
            Matrix matrix = svd.getRotation();
            Atom shift = svd.getTranslation();

            for(Atom a : ca2aligned) {
               Calc.rotate(a, matrix);
               Calc.shift(a, shift);
            }

            return SVDSuperimposer.getTMScore(ca1aligned, ca2aligned, ca1.length, ca2.length);
         } catch (StructureException e) {
            e.printStackTrace();
            return Double.NaN;
         }
      }

      @Override
      public String getName() {
         return "Ref_TM";
      }
   }

   public static class Alignment extends TMScoreMetric {

      @Override
      public double calculate(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {

         try {
            return AFPChainScorer.getTMScore(align,ca1,ca2);
         } catch (StructureException e){
            return Double.NaN;
         }
      }

      @Override
      public String getName() {
         return "Aln_TM";
      }
   }
   
   @Override
   public String format(double result) {
      return String.format("%.4f", result);
   }
}
