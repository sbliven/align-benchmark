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
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FastaStructureParser;
import org.biojava.bio.structure.io.StructureSequenceMatcher;
import org.biojava.bio.structure.scop.RemoteScopInstallation;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava.bio.structure.scop.ScopInstallation;
import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.Profile;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.MultipleSequenceAlignment;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.CasePreservingProteinSequenceCreator;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.biojava3.core.sequence.io.ProteinSequenceCreator;
import org.biojava3.core.sequence.io.template.FastaHeaderParserInterface;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;
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
	 * If true, consider upper-case letters aligned and lower-case letters
	 * unaligned.
	 */
	//private boolean caseSensitive = false;
	
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
		
		String filename;
		FastaHeaderParserInterface<ProteinSequence, AminoAcidCompound> headerParser;
		try {
			
			filename = "/home/sbliven/Documents/benchmark/Scheeff_kinase_align.fasta";
			
			// parse headers like '> d4hhba_'
			headerParser = new FastaHeaderParserInterface<ProteinSequence, AminoAcidCompound>() {
				public void parseHeader(String header,
						ProteinSequence sequence) {
					String scopDomain = header.trim().substring(0, 7);
					sequence.setOriginalHeader(header);
					sequence.setAccession(new AccessionID(scopDomain));
				}
			};

			filename = "/Users/blivens/dev/bourne/benchmarks/youkha/1KQ1-1C4Q.a2m";
			// >pdb|1KQ1|A|gi|22219069 [Staphylococcus aureus] 1.55 A Crystal Structure Of 
			//headerParser = new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>();
			headerParser = new FastaHeaderParserInterface<ProteinSequence, AminoAcidCompound>() {
				public void parseHeader(String header,
						ProteinSequence sequence) {
					sequence.setOriginalHeader(header);
					header = header.trim();
					String[] fields = header.split("\\|");
					
					// PDBID+CHAIN
					String accession = fields[1]+"."+fields[2];
					
					sequence.setAccession(new AccessionID(accession));
				}
			};

			
			AminoAcidCompoundSet aaSet;
			aaSet = AminoAcidCompoundSet.getAminoAcidCompoundSet();
			
			FastaReader<ProteinSequence, AminoAcidCompound> reader
				= new FastaReader<ProteinSequence, AminoAcidCompound>(
					new File(filename),
					headerParser,
					new CasePreservingProteinSequenceCreator(aaSet));
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
			// Requires using a local installation rather than RemoteScopInstallation
			ScopInstallation scop = new ScopInstallation();
			scop.setScopVersion("1.65");
			ScopFactory.setScopDatabase(scop);

			AtomCache cache = new AtomCache();
			cache.setFetchFileEvenIfObsolete(true);
			cache.setStrictSCOP(false);
			
			
			
			// order the sequences
			List<String> sequenceIDs = new ArrayList<String>(sequences.keySet());
			
			
			
			// Match each sequence  to a series of PDB Residue numbers
			List<List<PDBResidue>> residues = new ArrayList<List<PDBResidue>>();
			//List<List<Boolean>> isGap = new ArrayList<List<Boolean>>(); // Keep track of gaps
			for(String accession : sequenceIDs) {
				if(accession.equals("d1cjaa_")) break;

				ProteinSequence seq = sequences.get(accession);
								
				System.out.println("Fetching "+accession);
				Structure struct = cache.getStructure(accession);
				List<PDBResidue> match = PDBResidue.matchSequenceToStructure(seq, struct);
				
				assert( match.size() == seq.getLength());
				
				residues.add(match);
				
				// check for gaps
				Collection<Object> cases = seq.getUserCollection();// set by CaseSensitiveProteinSequenceCreator
				assert(cases.size() == match.size());
				//List<Boolean> gaps = new ArrayList<Boolean>(cases.size());
				int i=0;
				for( Object res : cases ) {
					boolean gap = res.equals(Character.UPPERCASE_LETTER);
					if(gap) {
						match.set(i, null);
					}
					//gaps.add( gap );
					i++;
				}
				//isGap.add(gaps);
			}
			
			// print the AA sequence and corresponding residues
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

	/**
	 * Removes all gaps from a protein sequence
	 * @param gapped
	 * @return
	 */
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
	/*
	/**
	 * @return whether to consider upper-case letters aligned and lower-case
	 *  letters unaligned.
	 *
	public boolean isCaseSensitive() {
		return caseSensitive;
	}

	/**
	 * @param caseSensitive If true, consider upper-case letters aligned and lower-case letters
	 * unaligned.
	 *
	public void setCaseSensitive(boolean caseSensitive) {
		this.caseSensitive = caseSensitive;
	}

	/*
	private static class CaseInsensitiveAminoAcidCompoundSet extends AminoAcidCompoundSet {
		public CaseInsensitiveAminoAcidCompoundSet() {
			super();
			
			List<AminoAcidCompound> upperCaseAAs = this.getAllCompounds();
			//List<AminoAcidCompound> lowerCaseAAs = new ArrayList<AminoAcidCompound>(upperCaseAAs.size());
			
			for(AminoAcidCompound aa:upperCaseAAs) {
				String upper = aa.getShortName();
				String lower = upper.toLowerCase();
				AminoAcidCompound aa2 = new AminoAcidCompound(this, lower, aa.getLongName(), aa.getDescription(), aa.getMolecularWeight());
				
				//lowerCaseAAs.add(aa2);
				aminoAcidCompoundCache.put(lower, aa2);
			}
		}
		
	    private final static CaseInsensitiveAminoAcidCompoundSet aminoAcidCompoundSet = new CaseInsensitiveAminoAcidCompoundSet();

	    public static CaseInsensitiveAminoAcidCompoundSet getAminoAcidCompoundSet() {
	        return aminoAcidCompoundSet;
	    }
	}
	*/
	
	public static void main2(String[] args)
	{
		
		String filename;
		FastaHeaderParserInterface<ProteinSequence, AminoAcidCompound> headerParser;
		try {
			
			filename = "/home/sbliven/Documents/benchmark/Scheeff_kinase_align.fasta";
			filename = "/Users/blivens/dev/bourne/benchmarks/Scheeff_kinase_align.fasta";
			// parse headers like '> d4hhba_'
			headerParser = new FastaHeaderParserInterface<ProteinSequence, AminoAcidCompound>() {
				public void parseHeader(String header,
						ProteinSequence sequence) {
					String scopDomain = header.trim().substring(0, 7);
					sequence.setOriginalHeader(header);
					sequence.setAccession(new AccessionID(scopDomain));
				}
			};

			
			filename = "/Users/blivens/dev/bourne/benchmarks/youkha/1KQ1-1C4Q.a2m";
			// >pdb|1KQ1|A|gi|22219069 [Staphylococcus aureus] 1.55 A Crystal Structure Of 
			//headerParser = new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>();
			headerParser = new FastaHeaderParserInterface<ProteinSequence, AminoAcidCompound>() {
				public void parseHeader(String header,
						ProteinSequence sequence) {
					sequence.setOriginalHeader(header);
					header = header.trim();
					String[] fields = header.split("\\|");
					
					// PDBID+CHAIN
					String accession = fields[1]+"."+fields[2];
					
					sequence.setAccession(new AccessionID(accession));
				}
			};
			
			
			
			// Use SCOP 1.65 (Dec 2003)
			// Requires using a local installation rather than RemoteScopInstallation
			ScopInstallation scop = new ScopInstallation();
			scop.setScopVersion("1.65");
			ScopFactory.setScopDatabase(scop);

			AtomCache cache = new AtomCache();
			cache.setFetchFileEvenIfObsolete(true);
			cache.setStrictSCOP(false);
			
			AminoAcidCompoundSet aaSet = AminoAcidCompoundSet.getAminoAcidCompoundSet();
			
			// parse file
			FastaStructureParser parser = new FastaStructureParser(
					new File(filename),
					headerParser,
					new CasePreservingProteinSequenceCreator(aaSet),
					cache);
			parser.process();
			
			ResidueNumber[][] residues = parser.getResidues();
			ProteinSequence[] sequences = parser.getSequences();
			Structure[] structures = parser.getStructures();
			String[] accessions = parser.getAccessions();
			
			// Set lowercase residues to null too
			for(int structNum = 0; structNum<sequences.length;structNum++) { // should iterate in the same order
				ProteinSequence seq = sequences[structNum];
				Collection<Object> userCollection = seq.getUserCollection();
				assert(userCollection != null); // should have been set by seq creator
				assert(userCollection.size() == residues[structNum].length);
				int pos = 0;
				for(Object isAligned : userCollection) {
					assert(isAligned instanceof Boolean);
					if(!(Boolean)isAligned) {
						residues[structNum][pos] = null;
					}
					pos++;
				}
			}

			
			// Remove alignment columns with a gap
			residues = StructureSequenceMatcher.removeGaps(residues);

			System.out.println();
			for(int i=0; i<residues.length;i++) {
				Structure struct = structures[i];
				System.out.print(accessions[i]+":\t");
				for(int j=0;j<residues[i].length;j++) {
					ResidueNumber res = residues[i][j];
					String name;
					if(res == null) {
						name = "-";
					} else {
						Group g = struct.getChainByPDB(res.getChainId()).getGroupByPDB(res);
						AminoAcidCompound aa = aaSet.getCompoundForString(g.getPDBName());
						name = aa.getShortName();
					}
					System.out.print(name);
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
	}

	
}