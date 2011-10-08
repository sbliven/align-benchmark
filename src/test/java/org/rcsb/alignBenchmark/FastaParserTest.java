/**
 * 
 */
package org.rcsb.alignBenchmark;

import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;

import org.biojava3.core.sequence.ProteinSequence;


/**
 * @author sbliven
 *
 */
public class FastaParserTest extends TestCase{

	public void testRemoveGapsPS() {
		String ungappedSeq = "ACDEFGHIKLMNPQRSTVWY";
		String gappedSeq;
		
		ProteinSequence ungapped = new ProteinSequence(ungappedSeq);
		ProteinSequence gapped;
		ProteinSequence degapped;
		
		gappedSeq = "ACDE-FGHIKL---MNPQRS--TVWY";
		gapped = new ProteinSequence(gappedSeq);
		degapped = FastaParser.removeGaps(gapped);
		assertEquals("Failed to remove middle gaps",
				ungappedSeq, degapped.getSequenceAsString());
		
		gappedSeq = "-ACDEFGHIKLMNPQRSTVWY-";
		gapped = new ProteinSequence(gappedSeq);
		degapped = FastaParser.removeGaps(gapped);
		assertEquals("Failed to remove end gaps",
				ungappedSeq, degapped.getSequenceAsString());
	}
	
	
	public void testRemoveGapsRes() {
		Integer[][] ungappedResNums = new Integer[][] { 
				new Integer[] { 1,2,3,4,5, },
				new Integer[] { 1,2,3,4,5, },
				new Integer[] { 1,2,3,4,5, },
		};
		
		Integer[][] gappedResNums;
		
		gappedResNums = new Integer[][] {
				new Integer[] {    0,1,2,   0,3,   0,null,4,5,null },
				new Integer[] { null,1,2,   0,3,null,null,4,5,   0 },
				new Integer[] {    0,1,2,null,3,null,null,4,5,   0 },
		};
		
		List<List<PDBResidue>> ungapped = makePDBResidueLists(ungappedResNums);
		List<List<PDBResidue>> gapped = makePDBResidueLists(gappedResNums);
		
		List<List<PDBResidue>> degapped = FastaParser.removeGaps(gapped);
		
		// deep equals
		assertEquals("Wrong size.",ungapped.size(),degapped.size());
		for(int prot = 0; prot< ungapped.size();prot++) {
			List<PDBResidue> ungappedProt = ungapped.get(prot);
			List<PDBResidue> degappedProt = degapped.get(prot);
			
			assertEquals("Wrong size of prot "+prot,ungappedProt.size(),degappedProt.size());
			
			for(int pos = 0; pos<ungappedProt.size();pos++) {
				PDBResidue ungappedRes = ungappedProt.get(pos);
				PDBResidue degappedRes = degappedProt.get(pos);
				
				assertEquals("Res at pos "+pos+" of prot "+prot+" differ.",
						ungappedRes, degappedRes);
			}
		}
		
	}


	private List<List<PDBResidue>> makePDBResidueLists(
			Integer[][] arr) {
		
		List<List<PDBResidue>> outer = new ArrayList<List<PDBResidue>>(arr.length);
		for(int i=0;i<arr.length;i++) {
			List<PDBResidue> inner = new ArrayList<PDBResidue>(arr[i].length);
			for(int j=0;j<arr[i].length;j++) {
				PDBResidue res = null;
				if(arr[i][j] != null) {
					res = new PDBResidue(arr[i][j].toString());
				}
				inner.add(res);
			}
			outer.add(inner);
		}
		return outer;
	}
}
