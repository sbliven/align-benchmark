/**
 * 
 */
package org.rcsb.alignBenchmark;

import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.SeqRes2AtomAligner;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SimpleSubstitutionMatrix;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;

/**
 * Small bean to hold information about a single residue in the PDB
 * if we don't want to create a full {@link org.biojava.bio.structure.Group Group} object.
 * @author Spencer Bliven
 */
public class PDBResidue extends ResidueNumber{

	private static final long serialVersionUID = 2708507169769496311L;
	
	private String aaName; //3-letter code
	
	/**
	 * @param aaName 3-letter String giving the amino acid
	 * @param residueNum The {@link org.biojava.bio.structure.Group#getPDBCode() residue code}
	 *  for this residue (residue number + insertion code)
	 */
	public PDBResidue(ResidueNumber residueNum, String aaName) {
		super(residueNum.getChainId(), residueNum.getSeqNum(), residueNum.getInsCode());
		this.aaName = aaName;
	}
	
	/**
	 * @param residueNum The {@link org.biojava.bio.structure.Group#getPDBCode() residue code}
	 *  for this residue (residue number + insertion code)
	 */
	public PDBResidue(ResidueNumber residueNum) {
		this(residueNum,null);
	}
	
	/**
	 * Creates a PDBResidue from a string giving the number and insertion code.
	 * @param residueNum The {@link org.biojava.bio.structure.Group#getPDBCode() residue code}
	 *  for this residue (residue number + insertion code)
	 * @see ResidueNumber#fromString(String)
	 */
	public PDBResidue(String residueNum) {
		this(residueNum, null);
	}
	
	/**
	 * @param residueNum The {@link org.biojava.bio.structure.Group#getPDBCode() residue code}
	 *  for this residue (residue number + insertion code)
	 * @param aaName 3-letter String giving the amino acid
	 */
	public PDBResidue(String residueNum, String aaName) {
		this(ResidueNumber.fromString(residueNum),aaName);
	}

	/**
	 * Creates a PDBResidue from a given Group.
	 * 
	 * The PDB number, insertion code, and 3-letter name of the Group will be
	 * stored in the new PDBResidue.
	 * 
	 * @param g A Group, usually representing one amino acid of a protein
	 */
	public PDBResidue(Group g) {
		this(g.getResidueNumber(),g.getPDBName());
	}

	/**
	 * @return the 3-letter amino acid code, or null if none is set
	 */
	public String getAaName() {
		return aaName;
	}

	/**
	 * @param aaName the 3-letter amino acid code
	 */
	public void setAaName(String aaName) {
		this.aaName = aaName;
	}

	
	/**
	 * Prints the amino acid, if given, followed by the residue number, insertion code, and chain.
	 * <p>
	 * Examples:<pre>
	 * LEU.5.A
	 * 100F.G
	 * </pre>
	 */
	public String toDescriptionString() {
		String str = "";
		if(aaName != null) {
			str+=aaName+".";
		}
		return String.format("%s%s.%s", str,super.toString(),this.getChainId() );
	}

	/**
	 * @return
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result + ((aaName == null) ? 0 : aaName.hashCode());
		return result;
	}

	/**
	 * Requires that the other object also be a PDBResidue and that all member fields be equal
	 * @param obj
	 * @return
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (!super.equals(obj)) {
			return false;
		}
		if (!(obj instanceof PDBResidue)) {
			return false;
		}
		PDBResidue other = (PDBResidue) obj;
		if (aaName == null) {
			if (other.aaName != null) {
				return false;
			}
		} else if (!aaName.equals(other.aaName)) {
			return false;
		}
		return true;
	}
	
	//TODO The following two static methods (getProteinSequenceForStructure and
	//     matchSequenceToStructure) are not tightly coupled with the
	//     PDBResidue class. Consider placing them elsewhere. If moved,
	//     PDBResidueTests should also be renamed.
	
	/**
	 * Generates a ProteinSequence corresponding to the sequence of struct,
	 * and maintains a mapping from the sequence back to the original groups.
	 * 
	 * Chains are appended to one another. 'X' is used for heteroatoms.
	 * 
	 * @param struct
	 * @param groupIndexPosition An empty map to put
	 *  (residue index in returned ProteinSequence) -> (g
	 * @return
	 */
	public static ProteinSequence getProteinSequenceForStructure(Structure struct,
			Map<Integer,Group> groupIndexPosition) {
		
		if( groupIndexPosition != null) {
			groupIndexPosition.clear();
		}
		
		StringBuilder seqStr = new StringBuilder();
		
		for(Chain chain : struct.getChains()) {
			List<Group> groups = chain.getAtomGroups();
			Map<Integer,Integer> chainIndexPosition = new HashMap<Integer, Integer>();
			int prevLen = seqStr.length();
			
			// get the sequence for this chain
			String chainSeq = SeqRes2AtomAligner.getFullAtomSequence(groups, chainIndexPosition);
			seqStr.append(chainSeq);
			
			// fix up the position to include previous chains, and map the value back to a Group
			for(Integer seqIndex : chainIndexPosition.keySet()) {
				Integer groupIndex = chainIndexPosition.get(seqIndex);
				groupIndexPosition.put(prevLen + seqIndex, groups.get(groupIndex));
			}
		}
		
		
		return new ProteinSequence(seqStr.toString());
		
	}
	/**
	 * Given a sequence and the corresponding Structure, get the PDBResidue
	 * for each residue in the sequence.
	 * 
	 * <p>Smith-Waterman alignment is used to match the sequences. Residues
	 * in the sequence but not the structure or mismatched between sequence
	 * and structure will have a null PDBResidue, while
	 * residues in the structure but not the sequence are ignored.
	 * @param seq The protein sequence. Should match the sequence of struct very
	 * 	closely.
	 * @param struct The corresponding protein structure
	 * @return A list of PDBResidues of the same length as seq, containing either
	 *  the corresponding PDB location or null.
	 */
	public static List<PDBResidue> matchSequenceToStructure(ProteinSequence seq, Structure struct) {
		
		//1. Create ProteinSequence for struct while remembering to which group each residue corresponds
		Map<Integer,Group> atomIndexPosition   = new HashMap<Integer, Group>();
		
		ProteinSequence structSeq = getProteinSequenceForStructure(struct,atomIndexPosition);
		
		//2. Run Smith-Waterman to get the alignment
		// Identity substitution matrix with +1 for match, -1 for mismatch
		SubstitutionMatrix<AminoAcidCompound> matrix = 
			new SimpleSubstitutionMatrix<AminoAcidCompound>(
					AminoAcidCompoundSet.getAminoAcidCompoundSet(),
					(short)1, (short)-1 );
		matrix = new SimpleSubstitutionMatrix<AminoAcidCompound>(
				AminoAcidCompoundSet.getAminoAcidCompoundSet(),
				new InputStreamReader(
						SimpleSubstitutionMatrix.class.getResourceAsStream("/blosum100.txt")),
				"blosum100");
        SequencePair<ProteinSequence, AminoAcidCompound> pair = 
        	Alignments.getPairwiseAlignment(seq, structSeq,
                PairwiseSequenceAlignerType.GLOBAL, new SimpleGapPenalty(), matrix);
		
        //System.out.print(pair.toString());
        
		//3. Convert the alignment back to PDBResidues
        AlignedSequence<ProteinSequence,AminoAcidCompound> alignedSeq = pair.getQuery();
        AlignedSequence<ProteinSequence,AminoAcidCompound> alignedStruct = pair.getTarget();
		
        System.out.println(pair.toString(80));
        
        assert(alignedSeq.getLength() == alignedStruct.getLength());
        
        System.out.format("%d/min{%d,%d}\n", pair.getNumIdenticals(),alignedSeq.getLength()-alignedSeq.getNumGaps(),
        		alignedStruct.getLength()-alignedStruct.getNumGaps());
        
        
        List<PDBResidue> correspondingResidues = new ArrayList<PDBResidue>(seq.getLength());
        
        for( int pos = alignedSeq.getStart().getPosition(); pos <= alignedSeq.getEnd().getPosition(); pos++ ) { // 1-indexed
        	//skip missing residues from sequence. These probably represent alignment errors
        	if(alignedSeq.isGap(pos)) {
        		int structIndex = alignedStruct.getSequenceIndexAt(pos)-1;
        		assert(structIndex > 0);//should be defined since seq gap
        		
        		Group g = atomIndexPosition.get(structIndex);
        		
        		System.err.format("Warning: residue %s %s in the Structure has no corresponding amino acid in the sequence.\n",
        				g.getChainId(),
        				g.getResidueNumber().toString());
        		continue;
        	}
        	
        	PDBResidue res = null;
        	
        	
        	if(! alignedStruct.isGap(pos) ) {
        		int structIndex = alignedStruct.getSequenceIndexAt(pos)-1;//1-indexed
        		Group g = atomIndexPosition.get(structIndex);
        		res = new PDBResidue(g);
        	}
        	correspondingResidues.add(res);
        }
        return correspondingResidues;
	}
}
