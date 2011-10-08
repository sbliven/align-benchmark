/**
 * 
 */
package org.rcsb.alignBenchmark;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava.bio.structure.scop.ScopInstallation;
import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.Profile;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.MultipleSequenceAlignment;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.biojava3.core.sequence.io.ProteinSequenceCreator;
import org.biojava3.core.sequence.io.template.FastaHeaderParserInterface;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;


/**
 * Parses multiple alignments stored in FASTA format
 * @author Spencer Bliven
 *
 */
public class FastaParser implements MultipleAlignmentParser
{
	// Everything required for a FastaReader
	private MultipleAlignment ma;
	
	/**
	 * 
	 * @param filename
	 * @param headerParser Parses the description line of each record,
	 * @param sequenceCreator
	 * @throws IOException 
	 */
	public FastaParser(String filename,
			FastaHeaderParserInterface<ProteinSequence, AminoAcidCompound> headerParser,
			SequenceCreatorInterface<AminoAcidCompound> sequenceCreator, boolean pairwiseMAs) throws IOException {
		this.ma = createMultipleAlignment(filename, headerParser, sequenceCreator, pairwiseMAs);
	}

	public FastaParser(String filename, boolean pairwiseMAs) throws IOException {
		this(filename,
				new SCOPFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
        		new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()), pairwiseMAs);
	}
	
	public FastaParser(String filename) throws IOException {
		this(filename,true);
	}
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args)
	{
		
		String filename = "/home/sbliven/Documents/benchmark/Scheeff_kinase_align.fasta";
		try {
			
			FastaHeaderParserInterface<ProteinSequence, AminoAcidCompound> headerParser;
			//headerParser = new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>();
			headerParser = new FastaHeaderParserInterface<ProteinSequence, AminoAcidCompound>() {
					public void parseHeader(String header,
							ProteinSequence sequence) {
						String scopDomain = header.trim().substring(0, 7);
						sequence.setOriginalHeader(header);
						sequence.setAccession(new AccessionID(scopDomain));
					}
				};

			FastaReader<ProteinSequence, AminoAcidCompound> reader
				= new FastaReader<ProteinSequence, AminoAcidCompound>(
					new File(filename),
					headerParser,
					new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
			LinkedHashMap<String,ProteinSequence> sequences = reader.process();

			//Profile<ProteinSequence, AminoAcidCompound> msa = new SimpleProfile<ProteinSequence,AminoAcidCompound>(null,null);
			
			/*			
			MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound> msa =
				new MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound>();
			int i = 0;
			for(ProteinSequence seq : sequences.values()) {
				msa.addAlignedSequence(seq);
				i++;
				if(i>=2) break;
			}
			for(ProteinSequence seq : msa.getAlignedSequences()) {
				System.out.format("Key=%s\nValue=%s\n",seq.getAccession(),seq.getSequenceAsString());
				AminoAcidCompound aa = seq.getCompoundAt(9);
				System.out.println(aa.getShortName());
			}
			*/
			
			// Use SCOP 1.65 (Dec 2003)
			ScopInstallation scop = new ScopInstallation();
			scop.setScopVersion("1.65");
			ScopFactory.setScopDatabase(scop);

			AtomCache cache = new AtomCache();
			cache.setFetchFileEvenIfObsolete(true);
			cache.setStrictSCOP(false);
			
			AminoAcidCompoundSet aaSet = AminoAcidCompoundSet.getAminoAcidCompoundSet();
			
			
			// order the sequences
			List<String> sequenceIDs = new ArrayList<String>(sequences.keySet());
			
			List<List<PDBResidue>> residues = new ArrayList<List<PDBResidue>>();
			for(String accession : sequenceIDs) {
				if(accession.equals("d1cjaa_")) break;

				ProteinSequence seq = sequences.get(accession);
				//String accession = seq.getAccession().toString();
				System.out.println("Fetchin "+accession);
				Structure struct = cache.getStructure(accession);
				List<PDBResidue> match = PDBResidue.matchSequenceToStructure(seq, struct);
				
				assert( match.size() == seq.getLength());
				
				residues.add(match);
				
			}
			
			for(int j = 0;j<residues.size();j++) {
				ProteinSequence seq = sequences.get(sequenceIDs.get(j));
				List<PDBResidue> match = residues.get(j);
				
				assert( match.size() == seq.getLength());
				
				for(int i = 1; i<=seq.getLength(); i++) { // 1-indexed
					System.out.format("%5s", seq.getCompoundAt(i).getShortName());
				}
				System.out.println();
				for(int i = 0; i<match.size();i++) {
					PDBResidue res = match.get(i);
					if(res != null) {
						String aaName = getShortName(res);
						System.out.format("%4s%1s", res.toString(),aaName);
					}
					else {
						System.out.format("%5s","-");
					}
				}
				System.out.println();
				
			}
			
			// Remove alignment columns with a gap
			residues = removeGaps(residues);

			System.out.println();
			for(List<PDBResidue> reslist : residues) {
				int i=0;
				for(PDBResidue res : reslist) {
//					System.out.format("%d\t%s",i,res);
//					System.out.println();
					System.out.print(getShortName(res));
					i++;
				}
				System.out.println();
			}
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (StructureException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		/*FastaParser parser = new FastaParser(filename);
		for(MultipleAlignment ma : parser) {
			System.out.println(ma.display());
		}*/
	}

	private static String getShortName(	PDBResidue res) {
		AminoAcidCompoundSet aaSet = AminoAcidCompoundSet.getAminoAcidCompoundSet();
		String aaName = "?";
		if(res.getAaName() != null) {
			AminoAcidCompound aa = aaSet.getCompoundForString(res.getAaName());
			aaName = aa.getShortName();
		}
		return aaName;
	}

	public static ProteinSequence removeGaps(ProteinSequence gapped) {
		final String gapString = "-";
		
		StringBuilder seq = new StringBuilder();
		
		AminoAcidCompound gap = gapped.getCompoundSet().getCompoundForString(gapString);
		
		for(int i=1; i<=gapped.getLength();i++) { //1-indexed
			AminoAcidCompound aa = gapped.getCompoundAt(i);
			if(! aa.equals(gap)) {
				seq.append(aa.getShortName());
			}
		}
		
		ProteinSequence ungapped = new ProteinSequence(seq.toString());
		
		return ungapped;
	}
	
	/**
	 * Creates a new list consisting of all columns of gapped where no row 
	 * contained a null value.
	 * 
	 * Here, "row" refers to the first index and "column" to the second, eg
	 * gapped.get(row).get(column)
	 * @param gapped A rectangular matrix containing null to mark gaps
	 * @return A new List without columns containing nulls
	 */
	public static List<List<PDBResidue>> removeGaps(List<List<PDBResidue>> gapped) {
		if(gapped == null) return null;
		
		final int nProts = gapped.size();
		List<List<PDBResidue>> ungapped = new ArrayList<List<PDBResidue>>(nProts);

		if(nProts < 1)
			return ungapped;
		
		
		int protLen = gapped.get(0).size();
		for(int i=0;i<nProts;i++) {
			int currProtLen = gapped.get(i).size();
			if(currProtLen != protLen) {
				throw new IllegalArgumentException(String.format(
						"Expected a square array, but row 0 has %d elements " +
						"while row %d has %d.", protLen,i,currProtLen));
				
			}
			ungapped.add(new ArrayList<PDBResidue>());
		}
		
		for(int pos=0;pos<protLen;pos++) {
			// Check if any rows have gaps
			boolean gap = false;
			for(List<PDBResidue> prot : gapped) {
				if(prot.get(pos) == null) {
					gap = true;
				}
			}
			
			if(!gap) {
				// Not a gap, so copy over all the PDBResidues
				for(int protNum = 0; protNum < nProts; protNum++) {
					ungapped.get(protNum).add(gapped.get(protNum).get(pos));
				}
			}
		}
		
		return ungapped;
	}
	/**
	 * Returns an iterator over the single multiple alignment represented by this Fasta file.
	 * @return
	 * @see java.lang.Iterable#iterator()
	 */
	public Iterator<MultipleAlignment> iterator()
	{
		List<MultipleAlignment> maList = Arrays.asList(ma);
		return maList.iterator();
	}

	/**
	 * Create a multiple alignment from a Fasta File
	 * @param filename
	 * @param headerParser
	 * @param sequenceCreator
	 * @param pairwiseMAs
	 * @return
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	private static MultipleAlignment createMultipleAlignment(String filename,
			FastaHeaderParserInterface<ProteinSequence, AminoAcidCompound> headerParser,
			SequenceCreatorInterface<AminoAcidCompound> sequenceCreator,
			boolean pairwiseMAs)
	throws FileNotFoundException, IOException
	{

		FastaReader<ProteinSequence, AminoAcidCompound> reader = new FastaReader<ProteinSequence, AminoAcidCompound>(
				new File(filename),
				headerParser,
				sequenceCreator);
		LinkedHashMap<String,ProteinSequence>sequences = reader.process();
		
		MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound> msa =
			new MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound>();
		for(ProteinSequence seq : sequences.values()) {
			msa.addAlignedSequence(seq);
		}
		
		return null; //TODO Stub
	}



	/**
	 * Converts scop domain identifiers (eg 'd1lnlb1') into a PDB ID and chain
	 * (eg '1lnl.B').
	 * @param scopID
	 * @return the extracted pdbid and chain, or null if the scopID is malformed.
	 *
	public static String getPDBName(String scopID) {
		Matcher match = scopRegex.matcher(scopID);
		if(!match.matches()) {
			return null;
		}
		if(!match.group(2).equals("_")) {
			return match.group(1)+"."+match.group(2).toUpperCase();
		} else {
			return match.group(1);
		}
	}
	/**
	 * Parses an array of scop domain identifiers, calling getPDBName for each.
	 * @param scopIDs
	 * @return
	 * @see getPDBName
	 *
	public static String[] getPDBNames(String[] scopIDs) {
		String[] pdbNames = new String[scopIDs.length];
		for(int i=0;i<scopIDs.length;i++) {
			pdbNames[i] = getPDBName(scopIDs[i]);
		}
		return pdbNames;
	}
	
	
	/**
	 * Parses a RIPC alignment file.
	 * <p>
	 * Format:<br/><pre>
	 * alignment := comment* labelPair resPair+ '\n'
	 * comment := '#' .* '\n' | '\n'
	 * labelPair := '#' (label1) '-' (label2)
	 * resPair := res ' ' res '\n'
	 * res := (3-letter amino acid code)'.'(pdb res number)'.'(insertion code)'.'(chain)</pre>
	 * Where outputs are listed in parentheses. Note that alignments must be separated by empty lines.
	 * <p>
	 * Example:<br/><pre>
	 * #d1hcy_2-d1lnlb1
	 * HIS.194._._	HIS.41._.B
	 * HIS.198._._	HIS.61._.B
	 * HIS.224._._	HIS.70._.B
	 * HIS.344._._	HIS.181._.B
	 *
	 * @author Spencer Bliven
	 *
	 *
	protected static class RIPCIterator implements Iterator<MultipleAlignment> {
		private BufferedReader ripc;
		private String[] nextLabels;

		private static final Pattern labelRegex =
			Pattern.compile("^\\#([^-]+)-([^-]+)$");
		private static final Pattern pairRegex = 
			Pattern.compile("^([A-Z]{3})\\.(\\d+)\\.(.)\\.(.)\\s+([A-Z]{3})\\.(\\d+)\\.(.)\\.(.)$");
		private static final Pattern commentRegex = 
			Pattern.compile("^(?:#.*|$)"); // labels are a subset of comments; check for labels first

		public RIPCIterator(String filename) throws IOException {
			this(new BufferedReader(new FileReader(filename)));
		}
		public RIPCIterator(BufferedReader ripc) throws IOException {
			this.ripc = ripc;
			nextLabels = null;
			skipComments();
		}
		//@Override
		public MultipleAlignment next()
		{
			List<List<PDBResidue>> residues = new ArrayList<List<PDBResidue>>(2);
			residues.add(new LinkedList<PDBResidue>());
			residues.add(new LinkedList<PDBResidue>());


			String line;
			try {
				line = ripc.readLine();

				while( line!=null ) {
					line = line.trim();
					// Check if line defines a residue pair
					Matcher pair = pairRegex.matcher(line);
					if(pair.matches()) {
//						String aa1 = pair.group(1);
//						String aa2 = pair.group(5);

						String pdb1 = pair.group(2);
						String pdb2 = pair.group(6);

						//Add insertion codes
						if( !pair.group(3).equals("_") ) {
							pdb1 += pair.group(3);
						}
						if( !pair.group(7).equals("_") ) {
							pdb2 += pair.group(7);
						}

//						String chain1 = pair.group(4);
//						String chain2 = pair.group(8);


						residues.get(0).add(new PDBResidue(pdb1));
						residues.get(1).add(new PDBResidue(pdb2));
					}
					else {
						// Check if line starts a new set of proteins
						Matcher labels=labelRegex.matcher(line);
						if(labels.matches()) {
							String[] names = nextLabels.clone();
							// To use full chains:
							//names = RIPCParser.getPDBNames(nextLabels);
							MultipleAlignment m = new MultipleAlignment(names,residues);

							nextLabels[0] = labels.group(1);
							nextLabels[1] = labels.group(2);

							return m;
						}
						else if(commentRegex.matcher(line).matches()) {
							// ignore comments
						}
						else {
							// Oops! Something's wrong
							throw new IllegalStateException("Formatting error. Unrecognized line:\n"+line);
						}
					}

					line = ripc.readLine();
				}
			} catch (IOException e) {
				throw new NoSuchElementException("IOException occured");
			}

			MultipleAlignment m = new MultipleAlignment(FastaParser.getPDBNames(nextLabels),residues);
			nextLabels = null;
			return m;
		}

		/**
		 * Skips the first level of comments and initializes nextLabels to the first alignment.
		 * @throws IOException 
		 * @throws IOException 
		 *
		private void skipComments() throws IOException {
			String line;
			line = ripc.readLine();

			while( line!=null ) {
				line = line.trim();
				Matcher labels=labelRegex.matcher(line);
				if(labels.matches()) {
					nextLabels = new String[2];
					nextLabels[0] = labels.group(1);
					nextLabels[1] = labels.group(2);
					return;
				}
				Matcher comments=commentRegex.matcher(line);
				if(!comments.matches()) {
					// Oops! Something's wrong
					throw new IllegalStateException("Formatting error. Expected comment or label, found:\n"+line);
				}

				line = ripc.readLine();
			}

			throw new NoSuchElementException();
		}
		

		//@Override
		public boolean hasNext()
		{
			return nextLabels != null;
		}



		/**
		 * Not implemented.
		 *
		//@Override
		public void remove()
		{
			throw new UnsupportedOperationException();
		}

	}

*/
}
