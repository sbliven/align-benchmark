/**
 * 
 */
package org.rcsb.alignBenchmark;


import java.util.ArrayList;
import java.util.Arrays;

import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.io.template.FastaHeaderParserInterface;

/**
 * Parses FastaHeaders containing SCOP IDs.
 * <p>
 * Example: <pre>d1bo1a_/1-292</pre>
 * <p>
 * This type of header is used in the Scheeff kinase alignment.
 * 
 * @author Spencer Bliven <sbliven@ucsd.edu>
 *
 */
public class SCOPFastaHeaderParser<S extends ProteinSequence, C extends AminoAcidCompound> implements FastaHeaderParserInterface<S, C> {

	
	public void parseHeader(String header, S sequence) {
		sequence.setOriginalHeader(header);
		
		String[] parts = header.split("/");
		String domain = parts[0].trim();
		String range = parts[1].trim();
		sequence.setAccession(new AccessionID(domain));
		ArrayList<String> notes = sequence.getNotesList();
		notes.add(range);
		//sequence.setNotesList(notes);
	}

}
