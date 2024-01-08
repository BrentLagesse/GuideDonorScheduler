import unittest
from unittest.mock import patch
import main

class TestDNAInverter(unittest.TestCase):

    def test_inverter(self):
        """Test converting DNA to reverse complement"""
        dna = 'AGACCTTGATTCGTGCTGTTTCTTCTCCTCAA'
        expected = 'TTGAGGAGAAGAAACAGCACGAATCAAGGTCT'
        
        actual = main.invert_dna(dna)  
        self.assertEqual(actual, expected)

class TestPAMFinder(unittest.TestCase):

    @patch('main.get_dna') 
    def test_find_pams(self, mock_get_dna):  
        """Test finding PAM sites in sequence""" 
        mock_dna = "AACCGGTTAACCGG"
        mock_get_dna.return_value = None, [mock_dna], None
        
        expected_pams = ['AACCGG', 'AACCGG']        
        actual_pams = [mock_dna[loc:loc+3] for loc in main.get_locations(mock_dna)[0]]
        
        self.assertEqual(actual_pams, expected_pams)

class TestMutations(unittest.TestCase):

    @patch('main.get_dna')
    def test_number_mutations(self, mock_get_dna):
        """Test correct number of mutations generated"""
        mock_dna = "AACCGGTTAACCGG"  
        mock_locs = [3] 
        mock_get_dna.return_value = None, [mock_dna], None
        
        mutations = main.get_all_mutations(mock_locs, [], mock_dna, mock_dna)
        self.assertEqual(len(mutations), 1) # placeholder, should check against actual config
          
if __name__ == '__main__':
    unittest.main()