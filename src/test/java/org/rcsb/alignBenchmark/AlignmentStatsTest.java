package org.rcsb.alignBenchmark;
/**
 * 
 */


import java.io.IOException;
import java.util.ArrayList;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;
import org.biojava.bio.structure.align.util.AtomCache;
import org.rcsb.alignBenchmark.MultipleAlignment;
import org.rcsb.alignBenchmark.metrics.*;

import junit.framework.TestCase;

/**
 * @author Spencer Bliven
 *
 */
public class AlignmentStatsTest extends TestCase {
	MultipleAlignment ref;
	Atom[][] structures;
	double tolerance;

	@Override
	public void setUp() throws StructureException, IOException {
		AtomCache cache = new AtomCache(System.getProperty("java.io.tmpdir"),true);


		ref = createMA1();

		assertEquals("Wrong number of proteins in ref alignment.",2,ref.getNames().length);
		String name1 = ref.getNames()[0];
		String name2 = ref.getNames()[1];

		structures = new Atom[2][];
		structures[0] = cache.getAtoms(name1);
		structures[1] = cache.getAtoms(name2);
		
		tolerance = 1.e-10;

	}

	public void testCEAlignment() throws StructureException {
		StructureAlignment aligner = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
		AFPChain alignment = aligner.align(structures[0],structures[1]);

		
		ArrayList<Metric> metrics = new ArrayList<Metric>();
		metrics.add(new ProteinLengthMetric(0));
		metrics.add(new ProteinLengthMetric(1));
		metrics.add(new MetaDataMetric("alignmentTime","time"));
		metrics.add(new AlignmentLengthMetric.Reference());
		metrics.add(new AlignmentLengthMetric.Alignment());
		metrics.add(new RMSDMetric.Reference());
		metrics.add(new RMSDMetric.Alignment());
		metrics.add(new TMScoreMetric.Reference());
		metrics.add(new TMScoreMetric.Alignment());
		metrics.add(new PercentCorrectMetric());
		metrics.add(new ConsistencyMetric() );
		metrics.add(new ConsistencyMetric(4) );
		metrics.add(new GapLengthMetric.Reference());
		// For fewer than 100 gap residues, the integer part gives number of gaps
		// and the fractional part gives total gap length
		metrics.add(new GapLengthMetric.Reference(1.01,.01));
		metrics.add(new GapLengthMetric.Alignment());
		metrics.add(new GapLengthMetric.Alignment(1.01,.01));
		
		metrics.add(new JiaCEScoreMetric.Reference());
		metrics.add(new JiaCEScoreMetric.Alignment());

		
		AlignmentStats stats = new AlignmentStats(metrics,ref, alignment, structures[0],structures[1], null);

		int i=0;
		assertEquals("Wrong Length for protein 1.", 78., stats.getResult(i++),tolerance);
		assertEquals("Wrong Length for protein 2.", 77., stats.getResult(i++),tolerance);
		i++;//assertFalse("No time assigned.",Double.isNaN(stats.getResult(i++)));
		assertEquals("Wrong Ref length.",72., stats.getResult(i++),tolerance);
		assertEquals("Wrong Aln length.",71., stats.getResult(i++),tolerance);
		
		// Note that the following numbers are not hand calculated
		assertEquals("Wrong Ref RMSD.",2.3844852594668304, stats.getResult(i++),tolerance);
		assertEquals("Wrong Aln RMSD.",4.153834616231908, stats.getResult(i++),tolerance);
		assertEquals("Wrong Ref TM-Score.",.6456601492360861, stats.getResult(i++),tolerance);
		assertEquals("Wrong Aln TM-Score.",.41129832046341386, stats.getResult(i++),tolerance);
		
		assertEquals("Wrong PercentCorrectMetric.",0., stats.getResult(i++),tolerance);
		assertEquals("Wrong 0-consistency.",0., stats.getResult(i++),tolerance);
		assertEquals("Wrong 4-consistency.",0., stats.getResult(i++),tolerance);

		assertEquals("Wrong number of gaps.",1., stats.getResult(i++),tolerance);
		assertEquals("Wrong gap penalty.", 1.03, stats.getResult(i++),tolerance);
		assertEquals("Wrong number of gaps.",3., stats.getResult(i++),tolerance);
		assertEquals("Wrong gap penalty.", 3.10, stats.getResult(i++),tolerance);

		assertEquals("Wrong Ref ce_score.",2.3844852594668304/72.*(1.+1./72.), stats.getResult(i++),tolerance);
		assertEquals("Wrong Aln ce_score.",4.153834616231908/71.*(1.+3./71.), stats.getResult(i++),tolerance);

	}
	
	public void testSWAlignment() throws StructureException {
		StructureAlignment aligner = StructureAlignmentFactory.getAlgorithm(SmithWaterman3Daligner.algorithmName);
		AFPChain alignment = aligner.align(structures[0],structures[1]);

		
		ArrayList<Metric> metrics = new ArrayList<Metric>();
		metrics.add(new ProteinLengthMetric(0));
		metrics.add(new ProteinLengthMetric(1));
//		metrics.add(new MetaDataMetric("alignmentTime","time"));
		metrics.add(new AlignmentLengthMetric.Reference());
		metrics.add(new AlignmentLengthMetric.Alignment());
		
		metrics.add(new RMSDMetric.Reference());
		metrics.add(new RMSDMetric.Alignment());
		metrics.add(new TMScoreMetric.Reference());
		metrics.add(new TMScoreMetric.Alignment());
		
		metrics.add(new PercentCorrectMetric());
		metrics.add(new ConsistencyMetric() );
		metrics.add(new ConsistencyMetric(4) );

		metrics.add(new GapLengthMetric.Reference());
		// For fewer than 100 gap residues, the integer part gives number of gaps
		// and the fractional part gives total gap length
		metrics.add(new GapLengthMetric.Reference(1.01,.01));
		metrics.add(new GapLengthMetric.Alignment());
		metrics.add(new GapLengthMetric.Alignment(1.01,.01));
		
		metrics.add(new JiaCEScoreMetric.Reference());
		metrics.add(new JiaCEScoreMetric.Alignment());

		
		AlignmentStats stats = new AlignmentStats(metrics,ref, alignment, structures[0],structures[1], null);
		
		int i=0;
		assertEquals("Wrong Length for protein 1.", 78., stats.getResult(i++),tolerance);
		assertEquals("Wrong Length for protein 2.", 77., stats.getResult(i++),tolerance);
		assertEquals("Wrong Ref length.",72., stats.getResult(i++),tolerance);
		assertEquals("Wrong Aln length.",35., stats.getResult(i++),tolerance);
		
		// Note that the following numbers are not hand calculated
		assertEquals("Wrong Ref RMSD.",2.3844852594668304, stats.getResult(i++),tolerance);
		assertEquals("Wrong Aln RMSD.",2.5333315417625233, stats.getResult(i++),tolerance);
		assertEquals("Wrong Ref TM-Score.",.6456601492360861, stats.getResult(i++),tolerance);
		assertEquals("Wrong Aln TM-Score.",.29790624754159745, stats.getResult(i++),tolerance);

		assertEquals("Wrong PercentCorrectMetric.",100.*35./72., stats.getResult(i++),tolerance);
		assertEquals("Wrong 0-consistency.",35./72., stats.getResult(i++),tolerance);
		assertEquals("Wrong 4-consistency.",35./72., stats.getResult(i++),tolerance);
		
		assertEquals("Wrong number of gaps.",1., stats.getResult(i++),tolerance);
		assertEquals("Wrong gap penalty.", 1.03, stats.getResult(i++),tolerance);
		assertEquals("Wrong number of gaps.",0., stats.getResult(i++),tolerance);
		assertEquals("Wrong gap penalty.", 0., stats.getResult(i++),tolerance);
		
		assertEquals("Wrong Ref ce_score.",2.3844852594668304/72.*(1.+1./72.), stats.getResult(i++),tolerance);
		assertEquals("Wrong Aln ce_score.",2.5333315417625233/35., stats.getResult(i++),tolerance);

	}
	
	public void testCECPAlignment() throws StructureException {
		StructureAlignment aligner = StructureAlignmentFactory.getAlgorithm(CeCPMain.algorithmName);
		AFPChain alignment = aligner.align(structures[0],structures[1]);

		
		ArrayList<Metric> metrics = new ArrayList<Metric>();
		metrics.add(new ProteinLengthMetric(0));
		metrics.add(new ProteinLengthMetric(1));
//		metrics.add(new MetaDataMetric("alignmentTime","time"));
		metrics.add(new AlignmentLengthMetric.Reference());
		metrics.add(new AlignmentLengthMetric.Alignment());
		
		metrics.add(new RMSDMetric.Reference());
		metrics.add(new RMSDMetric.Alignment());
		metrics.add(new TMScoreMetric.Reference());
		metrics.add(new TMScoreMetric.Alignment());
		
		metrics.add(new PercentCorrectMetric());
		metrics.add(new ConsistencyMetric() );
		metrics.add(new ConsistencyMetric(4) );

		metrics.add(new GapLengthMetric.Reference());
		// For fewer than 100 gap residues, the integer part gives number of gaps
		// and the fractional part gives total gap length
		metrics.add(new GapLengthMetric.Reference(1.01,.01));
		metrics.add(new GapLengthMetric.Alignment());
		metrics.add(new GapLengthMetric.Alignment(1.01,.01));
		
		metrics.add(new JiaCEScoreMetric.Reference());
		metrics.add(new JiaCEScoreMetric.Alignment());

		
		AlignmentStats stats = new AlignmentStats(metrics,ref, alignment, structures[0],structures[1], null);

		int i=0;
		assertEquals("Wrong Length for protein 1.", 78., stats.getResult(i++),tolerance);
		assertEquals("Wrong Length for protein 2.", 77., stats.getResult(i++),tolerance);
		assertEquals("Wrong Ref length.",72., stats.getResult(i++),tolerance);
		assertEquals("Wrong Aln length.",77., stats.getResult(i++),tolerance);
		
		// Note that the following numbers are not hand calculated
		assertEquals("Wrong Ref RMSD.",2.3844852594668304, stats.getResult(i++),tolerance);
		assertEquals("Wrong Aln RMSD.",2.7051866319000473, stats.getResult(i++),tolerance);
		assertEquals("Wrong Ref TM-Score.",.6456601492360861, stats.getResult(i++),tolerance);
		assertEquals("Wrong Aln TM-Score.",.6573435138027048, stats.getResult(i++),tolerance);

		assertEquals("Wrong PercentCorrectMetric.",100., stats.getResult(i++),tolerance);
		assertEquals("Wrong 0-consistency.",1., stats.getResult(i++),tolerance);
		assertEquals("Wrong 4-consistency.",1., stats.getResult(i++),tolerance);
		
		assertEquals("Wrong number of gaps.",1., stats.getResult(i++),tolerance);
		assertEquals("Wrong gap penalty.", 1.03, stats.getResult(i++),tolerance);
		assertEquals("Wrong number of gaps.",0., stats.getResult(i++),tolerance);
		assertEquals("Wrong gap penalty.", 0., stats.getResult(i++),tolerance);

		assertEquals("Wrong Ref ce_score.",2.3844852594668304/72.*(1.+1./72.), stats.getResult(i++),tolerance);
		assertEquals("Wrong Aln ce_score.",2.7051866319000473/77., stats.getResult(i++),tolerance);

	}
	/**
	 * Create multiple alignment for d1nkl__ vs d1qdma1 from RIPC.
	 * Note d1qdma1 is 1qdm(A:1S-104S)
	 * 
	 * Alignments:
	 * 
	 * RIPC:    42-76    2-38
	 *          3S-37S   67S-103S
	 * CECP:    1-39	 40-76     77
	 *          66S-104S 1S-37S    65S
	 * SW:      4-38
	 *          69S-103S
	 * CE:      1-37     38-49     54-60    61-71    74-77
	 *          1S-27S   65S-76S   77S-83S  88S-98S  99S-102S
	 * @return
	 */
	private static MultipleAlignment createMA1() {


		String[] names = new String[] {"d1nkla_", "d1qdma1"};

		PDBResidue[][] residues = new PDBResidue[][] {
				new PDBResidue[] {
						new PDBResidue("42"),
						new PDBResidue("43"),
						new PDBResidue("44"),
						new PDBResidue("45"),
						new PDBResidue("46"),
						new PDBResidue("47"),
						new PDBResidue("48"),
						new PDBResidue("49"),
						new PDBResidue("50"),
						new PDBResidue("51"),
						new PDBResidue("52"),
						new PDBResidue("53"),
						new PDBResidue("54"),
						new PDBResidue("55"),
						new PDBResidue("56"),
						new PDBResidue("57"),
						new PDBResidue("58"),
						new PDBResidue("59"),
						new PDBResidue("60"),
						new PDBResidue("61"),
						new PDBResidue("62"),
						new PDBResidue("63"),
						new PDBResidue("64"),
						new PDBResidue("65"),
						new PDBResidue("66"),
						new PDBResidue("67"),
						new PDBResidue("68"),
						new PDBResidue("69"),
						new PDBResidue("70"),
						new PDBResidue("71"),
						new PDBResidue("72"),
						new PDBResidue("73"),
						new PDBResidue("74"),
						new PDBResidue("75"),
						new PDBResidue("76"),
						// gap of 77,78,1
						new PDBResidue("2"),
						new PDBResidue("3"),
						new PDBResidue("4"),
						new PDBResidue("5"),
						new PDBResidue("6"),
						new PDBResidue("7"),
						new PDBResidue("8"),
						new PDBResidue("9"),
						new PDBResidue("10"),
						new PDBResidue("11"),
						new PDBResidue("12"),
						new PDBResidue("13"),
						new PDBResidue("14"),
						new PDBResidue("15"),
						new PDBResidue("16"),
						new PDBResidue("17"),
						new PDBResidue("18"),
						new PDBResidue("19"),
						new PDBResidue("20"),
						new PDBResidue("21"),
						new PDBResidue("22"),
						new PDBResidue("23"),
						new PDBResidue("24"),
						new PDBResidue("25"),
						new PDBResidue("26"),
						new PDBResidue("27"),
						new PDBResidue("28"),
						new PDBResidue("29"),
						new PDBResidue("30"),
						new PDBResidue("31"),
						new PDBResidue("32"),
						new PDBResidue("33"),
						new PDBResidue("34"),
						new PDBResidue("35"),
						new PDBResidue("36"),
						new PDBResidue("37"),
						new PDBResidue("38")
				},
				new PDBResidue[] {
						new PDBResidue("3S"),
						new PDBResidue("4S"),
						new PDBResidue("5S"),
						new PDBResidue("6S"),
						new PDBResidue("7S"),
						new PDBResidue("8S"),
						new PDBResidue("9S"),
						new PDBResidue("10S"),
						new PDBResidue("11S"),
						new PDBResidue("12S"),
						new PDBResidue("13S"),
						new PDBResidue("14S"),
						new PDBResidue("15S"),
						new PDBResidue("16S"),
						new PDBResidue("17S"),
						new PDBResidue("18S"),
						new PDBResidue("19S"),
						new PDBResidue("20S"),
						new PDBResidue("21S"),
						new PDBResidue("22S"),
						new PDBResidue("23S"),
						new PDBResidue("24S"),
						new PDBResidue("25S"),
						new PDBResidue("26S"),
						new PDBResidue("27S"),
						new PDBResidue("28S"),
						new PDBResidue("29S"),
						new PDBResidue("30S"),
						new PDBResidue("31S"),
						new PDBResidue("32S"),
						new PDBResidue("33S"),
						new PDBResidue("34S"),
						new PDBResidue("35S"),
						new PDBResidue("36S"),
						new PDBResidue("37S"),
						//gap of 38S, 65S, 66S
						new PDBResidue("67S"),
						new PDBResidue("68S"),
						new PDBResidue("69S"),
						new PDBResidue("70S"),
						new PDBResidue("71S"),
						new PDBResidue("72S"),
						new PDBResidue("73S"),
						new PDBResidue("74S"),
						new PDBResidue("75S"),
						new PDBResidue("76S"),
						new PDBResidue("77S"),
						new PDBResidue("78S"),
						new PDBResidue("79S"),
						new PDBResidue("80S"),
						new PDBResidue("81S"),
						new PDBResidue("82S"),
						new PDBResidue("83S"),
						new PDBResidue("84S"),
						new PDBResidue("85S"),
						new PDBResidue("86S"),
						new PDBResidue("87S"),
						new PDBResidue("88S"),
						new PDBResidue("89S"),
						new PDBResidue("90S"),
						new PDBResidue("91S"),
						new PDBResidue("92S"),
						new PDBResidue("93S"),
						new PDBResidue("94S"),
						new PDBResidue("95S"),
						new PDBResidue("96S"),
						new PDBResidue("97S"),
						new PDBResidue("98S"),
						new PDBResidue("99S"),
						new PDBResidue("100S"),
						new PDBResidue("101S"),
						new PDBResidue("102S"),
						new PDBResidue("103S")
				}
		};
		
		MultipleAlignment align = new MultipleAlignment(names, residues);
		return align;
	}
}
