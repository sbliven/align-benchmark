package org.rcsb.alignBenchmark;
/**
 * 
 */


import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.util.AtomCache;
import org.rcsb.alignBenchmark.MultipleAlignment;

import junit.framework.TestCase;

/**
 * @author Spencer Bliven
 *
 */
public class MultipleAlignmentTest extends TestCase {
	private AtomCache cache;
	
	public MultipleAlignmentTest() {
		cache = null;
	}
	public void setUp() {
		cache = new AtomCache();
	}

	public void testGetAlignmentMatrix() {
		/* Test on the following alignment:
		 #d1hcy_2-d1lnlb1
		 HIS.194._._	HIS.41._.B
		 HIS.198._._	HIS.61._.B
		 HIS.224._._	HIS.70._.B
		 HIS.344._._	HIS.181._.B
		 */
		String[] pdbIDs = new String[] {"1hcy.B", "1lnl.B"};

		PDBResidue[][] residues = new PDBResidue[][] {
				new PDBResidue[] {
						new PDBResidue("194"), 
						new PDBResidue("198"), 
						new PDBResidue("224"),
						new PDBResidue("344")
				},
				new PDBResidue[] {
						new PDBResidue("41"), 
						new PDBResidue("61"), 
						new PDBResidue("70"),
						new PDBResidue("181")
				},
		};
		//Expected output
		int[][] resIndices = new int[][] {
				new int[] { 193, 197, 223, 343, },
				new int[] { 43, 63, 72, 183 },
		};


		try {
			MultipleAlignment align = new MultipleAlignment(pdbIDs, residues);

			List<Atom[]> structures = new ArrayList<Atom[]>();
			for(String pdb :pdbIDs) {
				structures.add(cache.getAtoms(pdb));
			}

			int[][] alignMat = align.getAlignmentMatrix(pdbIDs, structures);

			assertTrue(Arrays.deepEquals(alignMat,resIndices));

			//			if(! Arrays.deepEquals(alignMat,resIndices)) {
			//				String msg = "Error in getAlignmentMatrix(). Expected:\n";
			//				msg += Arrays.deepToString(resIndices);
			//				msg += "\nFound:;";
			//				msg += Arrays.deepToString(alignMat);
			//
			//				System.err.println(msg);
			//			}
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());
		}

	}

	public void testFindAlignment() {
		try {

			//Atoms for 1nls
			Atom[] struct = cache.getAtoms("1nls");

			//Some residues which should be found in 1nls, and their corresponding numbers
			HashMap<PDBResidue,Integer> residues = new HashMap<PDBResidue, Integer>();
			residues.put(new PDBResidue("8"), 7);
			residues.put(new PDBResidue("8","GLU"), 7);
			residues.put(new PDBResidue("8","ASP"), -1); // Ignores the amino acid

			//Use reflection to get around private function
			//calls findGroup()h
			Method findGroup = MultipleAlignment.class.getDeclaredMethod("findGroup",
					Atom[].class, PDBResidue.class, int.class);
			findGroup.setAccessible(true);

			for(PDBResidue res : residues.keySet()) {
				Integer result = (Integer) findGroup.invoke(null, struct, res, 0);
				assertEquals(String.format("%s:", res),residues.get(res), result);
			}
		} catch( Exception e ) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
}
