/**
 * 
 */
package org.rcsb.alignBenchmark;

import org.biojava.bio.structure.ResidueNumber;

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
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		String str = "";
		if(aaName != null) {
			str+=aaName+".";
		}
		return String.format("%s%s.%s", str,super.toString(),this.getChainId() );
	}
	
	
}
