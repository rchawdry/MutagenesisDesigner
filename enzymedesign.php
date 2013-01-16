<?php
    header('Content-Type: application/json');

    $dna = isset($_REQUEST['dna']) ? $_REQUEST['dna']:'';
    $mutationPos = isset($_REQUEST['mutationPos']) ? $_REQUEST['mutationPos']:'';
    $mutationLen = isset($_REQUEST['mutationLen']) ? $_REQUEST['mutationLen']:'';    
    $returnArray = isset($_REQUEST['returnArray']) ? $_REQUEST['returnArray']:'';    
    $maxCuts = 3; //Only return enzymes that cut at most 3 locations
    $minMatchSize = 4;
    $enzymeType = 0;
    $includeTypeIIs = 0;
    $includeTypeIIb = 0;
    
    //$dna='atggatttgatggacatggatgccatcaatgacctccacttctctggacacagctccctcaccaacacccctacttccgacgaggattatggttgcacctggaaccactggtctcccatcgtcaactgggacacctttacaggtgccccagatgactttcaccacctcatggacaccatcatcgaggacagaacaactgtcctcgaacagctcagtccatctatcaccacaaccaccaccaccaccacaaccacggatgaggaagaagaagaaatggaaacaacaaccacaacaacaacaacggcaataaaaacgcatgaagttggtgatgattcaaaagggctaaaactagtgcacctattgatggctggcgcggaagccttaaccggttcaaccaaaaaccgtgacttagctcgagtgatattgatccggctcaaggaattagtctcacaacatgctaacggctccaacatggagagacttgccgctcatttcaccgaagcccttcatggactactagaaggtgctgggggtgcgcataataaccaccatcaccacaacaacaacaaacattacctcaccacaaatggtccccacgacaaccaaaacgacacacttgctgctttccaactacttcaagacatgtcaccttacgtcaaattcggccacttcacagccaatcaagccatcatcgaagccgtggcccacgaacgccgcgtccatgtcatcgactacgatatcatggaaggggtccaatgggcctcattaatccaatcacttgcttccaacaacaacggtccacatcttcggataaccgcattatcgcggacaggaaccggccggaggtcaatcgccaccgtacaagaaaccggcaggcgattaacctcctttgccgcttccctcggacaaccattttcttttcatcactgcagattagactccgatgaaactttccgtccttctgcattaaagcttgtacgtggagaggctttggttttcaactgcatgctgaatttacctcaccttagttaccgtgcaccggaatcagttgcttcgtttttaaacggagcaaaaacgttgaatccaaagcttgtgacgttggttgaagaagaagttggttctgttattggtgggtttgtggaaaggttcatggactcacttcatcattattcagcggtttttgattcattggaggctggttttccgatgcagaaccgggcccggactttggtggagagggtattttttgggcctagaattgcaggctcgttgggccggatatatagaacgggtggtgaagaggagagaaggtcatggggggagtggttaggtgaggtagggtttaggggagttccggtgagctttgcaaatcattgtcaagcaaagctgttgttagggctttttaatgatgggtatagggttgaagaggtgggtgtgggtagtaataagttggtgttggattggaaatcaagacgtttgctttctgcttccctttggacttgttcttcttcagattctgatttatag';
    if (strlen($dna)>0)
    {
        $e = new RestrictionEnzyme();
        
        $enzyme_list = $e->getEnzymeList($dna, $mutationPos, $mutationLen, $maxCuts);
        if (count($enzyme_list)>0)
        { 
            if ($returnArray)
                echo json_encode(array_values($enzyme_list));
            else
                echo json_encode($enzyme_list);
        }
    }
    else
    {
        echo "Please submit a DNA Sequence using the 'dna' form parameter";
    }

class RestrictionEnzyme {

    public function getEnzymeList($dna, $mutationPos, $mutationLen, $maxCuts, $minMatchSize = 4, $enzymeType = 0, $includeTypeIIs = 0, $includeTypeIIb = 0) {
        /*
         *  Gets a list of possible restriction enzymes given a given dna fragment
         *  $dna = string of DNA (ex GATACGA)
         *  $minMatchSize = minimum size for returned enzumes
         *  $enzymeType  = 0 - All enzymes (DEFAULT)
         *                 1 - Only Blunt-end enzymes
         *                 2 - Only Overhang enzymes
         */
      
        
        $sequence = preg_replace("/\W|\d/", "", strtoupper($dna));
        
        $enzymeLibrary = $this->getEnzymesTypeII();
        if ($includeTypeIIs) $enzymeLibrary = asort(array_merge($enzymeLibrary, $this->getEnzymesTypeIIs()));
        if ($includeTypeIIb) $enzymeLibrary = asort(array_merge($enzymeLibrary, $this->getEnzymesTypeIIb()));

        $enzymeLibrary = $this->removeEnzymes($enzymeLibrary, $minMatchSize, $enzymeType);
        $digestion=$this->digestDna($sequence,$enzymeLibrary);
        $enzymeList = $this->getEnzymeDetails($digestion, $enzymeLibrary);

        $this->assignRankings($enzymeList, $mutationPos, $mutationLen, strlen($dna));
        
        if ($maxCuts>0){
            foreach ($enzymeList as $key=>$value){
                if (count($value['cuts'])>$maxCuts){
                    unset($enzymeList[$key]);
                }
            }
        }
        
        uasort($enzymeList, function($a, $b){
            return $a['sort']<$b['sort'] ? -1 : 1;
            //return $a['rank']>$b['rank'] ? -1 : 1;
        });
        
        return $enzymeList;
        
        
    }
    
    private function getEnzymeDetails($selectedEnzymes, $enzymeLibrary){
        $enzymeList = array();
        foreach($selectedEnzymes as $enzyme=>$enzymeData){
            $enzymeList[$enzyme]['id']=  $enzyme;
            $enzymeList[$enzyme]['isomers']=  str_replace(',', ', ',$enzymeLibrary[$enzyme][0]);
            $enzymeList[$enzyme]['type']=($enzymeLibrary[$enzyme][5]==0 ? 'B':'O');
            $enzymeList[$enzyme]['pattern']=$enzymeLibrary[$enzyme][1];
            $enzymeList[$enzyme]['patternLength']=$enzymeLibrary[$enzyme][3];
            $enzymeList[$enzyme]['numCuts']=count($enzymeData['cuts']);
            $enzymeList[$enzyme]['cuts']=$enzymeData['cuts'];
            $enzymeList[$enzyme]['sort']=$enzymeData['sort'];
            $enzymeList[$enzyme]['cutList']=implode(', ',array_keys($enzymeData['cuts']));
        }
        return $enzymeList;
    }
    
    private function digestDna($sequence, $enzymes_array){
                       
        $digestion = array();
        foreach ($enzymes_array as $enzyme => $val){        
            
            // this is to put together results for IIb endonucleases, which are computed as "enzyme_name" and "enzyme_name@"
            $enzyme2=str_replace("@","",$enzyme);

            // split sequence based on pattern from restriction enzyme
            $fragments = preg_split("/".$enzymes_array[$enzyme][2]."/", $sequence,-1,PREG_SPLIT_DELIM_CAPTURE);
            reset ($fragments);
            $maxfragments=sizeof($fragments);

            // when sequence is cleaved ($maxfragments>1) start further calculations
            if ($maxfragments>1){
                $recognitionposition = strlen($fragments[0]);
                $counter_cleavages = 0;
                $list_of_cleavages = "";

                // for each frament generated, calculate cleavage position,
                //    add it to a list, and add 1 to counter
                $minCleavagePosition = -1;
                for ($i = 2; $i < $maxfragments; $i+=2) {
                    $cleavageposition = $recognitionposition + $enzymes_array[$enzyme][4];
                    
                    if ($cleavageposition<$minCleavagePosition || $minCleavagePosition==-1){
                        $minCleavagePosition = $cleavageposition;
                    }
                    $digestion[$enzyme2]["cuts"][$cleavageposition] = $cleavageposition;

                    /*
                    // As overlapping may occur for many endonucleases,
                    //   a subsequence starting in position 2 of fragment is calculate
                    $subsequence = substr($fragments[$i - 1], 1) . $fragments[$i] . substr($fragments[$i + 1], 0, 40);
                    $subsequence = substr($subsequence, 0, 2 * $enzymes_array[$enzyme][3] - 2);

                    //Previous process is repeated
                    // split subsequence based on pattern from restriction enzyme
                    $fragments_subsequence = preg_split($enzymes_array[$enzyme][2], $subsequence);

                    // when subsequence is cleaved start further calculations
                    if (sizeof($fragments_subsequence) > 1) {

                        // for each fragment of subsequence, calculate overlapping cleavage position,
                        //    add it to a list, and add 1 to counter
                        $overlapped_cleavage = $recognitionposition + 1 + strlen($fragments_subsequence[0]) + $enzymes_array[$enzyme][4];
                        $digestion[$enzyme2]["cuts"][$overlapped_cleavage] = "";
                    }
*/
                    // this is a counter for position
                    $recognitionposition+=strlen($fragments[$i - 1]) + strlen($fragments[$i]);
                }
                $digestion[$enzyme2]['sort']=$minCleavagePosition;
                //$digestion[$enzyme2]['cutList']=implode(',',asort(array_keys($digestion[$enzyme2]['cuts'])));
            }
        }
        return $digestion;
    }
    
    private function assignRankings(&$enzymes, $mutationPos, $mutationLen, $dnaLen){
        
        $minimumDistanceFromMutation = 150;
        
        $bestStartPos = ($mutationPos - $minimumDistanceFromMutation<0)?0:$mutationPos - $minimumDistanceFromMutation;        
        $bestEndPos = ($mutationPos + $mutationLen + $minimumDistanceFromMutation>$dnaLen) ? $dnaLen : $mutationPos + $mutationLen + $minimumDistanceFromMutation;
        
        foreach($enzymes as $enzyme=>&$enzymeData){

            $score_cuts = 1/$enzymeData['numCuts'];
            
            $score_distance=0;
            foreach($enzymeData['cuts'] as $cutLabel=>$cutPos){
                if ($cutPos<=$bestStartPos){
                    $distance = $bestStartPos - $cutPos;
                    $score_distance += 1-($distance/$bestStartPos);
                }
               elseif ($cutPos>=$bestEndPos){
                    $distance = $cutPos - $bestEndPos;
                    $score_distance += 1-($distance/$dnaLen);
               }               
            }
            $score_distance = $score_distance/$enzymeData['numCuts'];
            
            $score_type = ($enzymeData['type']=="O") ? 1:0;
               
             
            
            $rank = ($score_cuts + $score_distance + $score_type)/3;
            $enzymes[$enzyme]['rank'] = $rank;
        }
    }
    
    private function removeEnzymes($enzymes, $minMatchSize, $enzymeType) {

        $enzymes_array2 = array();
        
        foreach ($enzymes as $enzyme => $val) {

            // if retype==1 -> only Blund ends (continue for rest)
            if ($enzymeType == 1 and $enzymes[$enzyme][5] != 0)
                continue;

            // if retype==2 -> only Overhang end (continue for rest)
            if ($enzymeType == 2 and $enzymes[$enzyme][5] == 0)
                continue;

            // Only endonucleases with which recognized in template a minimum of bases (continue for rest)
            if ($minMatchSize > $enzymes[$enzyme][6])
                continue;

            $enzymes_array2[$enzyme] = $val;
        }

        return $enzymes_array2;
    }

    private function getEnzymesTypeII() {
        // this array includes all endonucleases related information required in this script
        // All enzymes with the same recognition pattern are grouped
        // The following information is provided for each endonuclease
        //  "AasI" => array(                             First endonuclease of the list
        //         0 => "AasI,DrdI,DseDI",               All endonucleases recognizing the same pattern
        //         1 => "GACNN_NN'NNGTC",                Recognition pattern
        //         2 => "(GAC......GTC)",                Recognition pattern for computing
        //         3 => 12,                              Length of all recognition pattern
        //         4 => 7,                               Cleavage position in upper strand
        //         5 => -2,                              Cleavage position in lower strand, relative to previous one
        //         6 => 6                                Number of non-N bases within recognition pattern
        //        ),

        $enzymes_array = Array(
            "AasI" => array("AasI,DrdI,DseDI", "GACNN_NN'NNGTC", "(GAC......GTC)", 12, 7, -2, 6),
            "AatI" => array("AatI,Eco147I,PceI,SseBI,StuI", "AGG'CCT", "(AGGCCT)", 6, 3, 0, 6),
            "AatII" => array("AatII", "G_ACGT'C", "(GACGTC)", 6, 5, -4, 6),
            "AbsI" => array("AbsI", "CC'TCGA_GG", "(CCTCGAGG)", 8, 2, 4, 8),
            "Acc16I" => array("Acc16I,AviII,FspI,NsbI", "TGC'GCA", "(TGCGCA)", 6, 3, 0, 6),
            "Acc65I" => array("Acc65I,Asp718I", "G'GTAC_C", "(GGTACC)", 6, 1, 4, 6),
            "AccB1I" => array("AccB1I,BanI,BshNI,BspT107I", "G'GYRC_C", "(GGCACC|GGCGCC|GGTACC|GGTGCC)", 6, 1, 4, 6),
            "AccB7I" => array("AccB7I,BasI,PflMI,Van91I", "CCAN_NNN'NTGG", "(CCA.....TGG)", 11, 7, -3, 6),
            "AccBSI" => array("AccBSI,BsrBI,MbiI", "CCG'CTC", "(CCGCTC|GAGCGG)", 6, 3, 0, 6),
            "AccI" => array("AccI,FblI,XmiI", "GT'MK_AC", "(GTCTAC|GTCGAC|GTATAC|GTAGAC)", 6, 2, 2, 6),
            "AccII" => array("AccII,Bsh1236I,BspFNI,BstFNI,BstUI,MvnI", "CG'CG", "(CGCG)", 4, 2, 0, 4),
            "AccIII" => array("AccIII,Aor13HI,BseAI,Bsp13I,BspEI,Kpn2I,MroI", "T'CCGG_A", "(TCCGGA)", 6, 1, 4, 6),
            "AciI" => array("AciI,BspACI,SsiI", "C'CG_C or G'CG_G", "(CCGC|GCGG)", 4, 1, 2, 4),
            "AclI" => array("AclI,Psp1406I", "AA'CG_TT", "(AACGTT)", 6, 2, 2, 6),
            "AcoI" => array("AcoI,CfrI,EaeI", "Y'GGCC_R", "(CGGCCA|CGGCCG|TGGCCA|TGGCCG)", 6, 1, 4, 6),
            "AcsI" => array("AcsI,ApoI,XapI", "R'AATT_Y", "(AAATTC|AAATTT|GAATTC|GAATTT)", 6, 1, 4, 6),
            "AcvI" => array("AcvI,BbrPI,Eco72I,PmaCI,PmlI,PspCI", "CAC'GTG", "(CACGTG)", 6, 3, 0, 6),
            "AcyI" => array("AcyI,BsaHI,BssNI,BstACI,Hin1I,Hsp92I", "GR'CG_YC", "(GACGCC|GACGTC|GGCGCC|GGCGTC)", 6, 2, 2, 6),
            "AdeI" => array("AdeI,DraIII", "CAC_NNN'GTG", "(CAC...GTG)", 9, 6, -3, 6),
            "AfaI" => array("AfaI,RsaI", "GT'AC", "(GTAC)", 4, 2, 0, 4),
            "AfeI" => array("AfeI,Aor51HI,Eco47III", "AGC'GCT", "(AGCGCT)", 6, 3, 0, 6),
            "AfiI" => array("AfiI,Bsc4I,BseLI,BsiYI,BslI", "CCNN_NNN'NNGG", "(CC.......GG)", 11, 7, -3, 4),
            "AflII" => array("AflII,BfrI,BspTI,Bst98I,MspCI,Vha464I", "C'TTAA_G", "(CTTAAG)", 6, 1, 4, 6),
            "AflIII" => array("AflIII", "A'CRYG_T", "(ACACGT|ACATGT|ACGCGT|ACGTGT)", 6, 1, 4, 6),
            "AgeI" => array("AgeI,AsiGI,BshTI,CspAI,PinAI", "A'CCGG_T", "(ACCGGT)", 6, 1, 4, 6),
            "AhdI" => array("AhdI,AspEI,BmeRI,DriI,Eam1105I,EclHKI", "GACNN_N'NNGTC", "(GAC.....GTC)", 11, 6, -1, 6),
            "AhlI" => array("AhlI,BcuI,SpeI", "A'CTAG_T", "(ACTAGT)", 6, 1, 4, 6),
            "AjiI" => array("AjiI,BmgBI,BtrI", "CAC'GTC", "(CACGTC|GACGTG)", 6, 3, 0, 6),
            "AjnI" => array("AjnI,EcoRII,Psp6I,PspGI", "'CCWGG_", "(CCAGG|CCTGG)", 5, 0, 5, 5),
            "AleI" => array("AleI,OliI", "CACNN'NNGTG", "(CAC....GTG)", 10, 5, 0, 6),
            "AluI" => array("AluI", "AG'CT", "(AGCT)", 4, 2, 0, 4),
            "Alw21I" => array("Alw21I,BsiHKAI,Bbv12I", "G_WGCW'C", "(GAGCAC|GAGCTC|GTGCAC|GTGCTC)", 6, 5, -4, 6),
            "Alw44I" => array("Alw44I,ApaLI,VneI", "G'TGCA_C", "(GTGCAC)", 6, 1, 4, 6),
            "AlwNI" => array("AlwNI,CaiI", "CAG_NNN'CTG", "(CAG...CTG)", 9, 6, -3, 6),
            "Ama87I" => array("Ama87I,AvaI,BmeT110I,BsiHKCI,BsoBI,Eco88I", "C'YCGR_G", "(CCCGAG|CCCGGG|CTCGAG|CTCGGG)", 6, 1, 4, 6),
            "ApaI" => array("ApaI", "G_GGCC'C", "(GGGCCC)", 6, 5, -4, 6),
            "ApeKI" => array("ApeKI,TseI", "G'CWG_C", "(GCAGC|GCTGC)", 5, 1, 3, 5),
            "AscI" => array("AscI,PalAI,SgsI", "GG'CGCG_CC", "(GGCGCGCC)", 8, 2, 4, 8),
            "AseI" => array("AseI,PshBI,VspI", "AT'TA_AT", "(ATTAAT)", 6, 2, 2, 6),
            "AsiSI" => array("AsiSI,RgaI,SgfI", "GCG_AT'CGC", "(GCGATCGC)", 8, 5, -2, 8),
            "Asp700I" => array("Asp700I,MroXI,PdmI,XmnI", "GAANN'NNTTC", "(GAA....TTC)", 10, 5, 0, 6),
            "AspA2I" => array("AspA2I,AvrII,BlnI,XmaJI", "C'CTAG_G", "(CCTAGG)", 6, 1, 4, 6),
            "AspI" => array("AspI,PflFI,PsyI,Tth111I", "GACN'N_NGTC", "(GAC...GTC)", 9, 4, 1, 6),
            "AspLEI" => array("AspLEI,BstHHI,CfoI,HhaI", "G_CG'C", "(GCGC)", 4, 3, -2, 4),
            "AspS9I" => array("AspS9I,BmgT120I,Cfr13I,Sau96I", "G'GNC_C", "(GG.CC)", 5, 1, 3, 4),
            "AssI" => array("AssI,BmcAI,ScaI,ZrmI", "AGT'ACT", "(AGTACT)", 6, 3, 0, 6),
            "AsuC2I" => array("AsuC2I,BcnI,BpuMI,NciI", "CC'S_GG", "(CCGGG|CCCGG)", 5, 2, 1, 5),
            "AsuNHI" => array("AsuNHI,BspOI,NheI", "G'CTAG_C", "(GCTAGC)", 6, 1, 4, 6),
            "AvaII" => array("AvaII,Bme18I,Eco47I,SinI,VpaK11BI", "G'GWC_C", "(GGACC|GGTCC)", 5, 1, 3, 5),
            "AxyI" => array("AxyI,Bse21I,Bsu36I,Eco81I", "CC'TNA_GG", "(CCT.AGG)", 7, 2, 3, 6),
            "BalI" => array("BalI,MlsI,MluNI,MscI,Msp20I", "TGG'CCA", "(TGGCCA)", 6, 3, 0, 6),
            "BamHI" => array("BamHI", "G'GATC_C", "(GGATCC)", 6, 1, 4, 6),
            "BanII" => array("BanII,Eco24I,EcoT38I,FriOI", "G_RGCY'C", "(GAGCCC|GAGCTC|GGGCCC|GGGCTC)", 6, 5, -4, 6),
            "BanIII" => array("BanIII,Bsa29I,BseCI,BshVI,BspDI,BspXI,Bsu15I,BsuTUI,ClaI", "AT'CG_AT", "(ATCGAT)", 6, 2, 2, 6),
            "BauI" => array("BauI", "C'ACGA_G'", "(CACGAG)", 6, 1, 4, 6),
            "BbeI" => array("BbeI", "G_GCGC'C", "(GGCGCC)", 6, 5, -4, 6),
            "BbuI" => array("BbuI,PaeI,SphI", "G_CATG'C", "(GCATGC)", 6, 5, -4, 6),
            "BbvCI" => array("BbvCI", "CC'TCA_GC or  GC'TGA_GG", "(CCTCAGC|GCTGAGG)", 7, 2, 3, 7),
            "BclI" => array("BclI,FbaI,Ksp22I", "T'GATC_A", "(TGATCA)", 6, 1, 4, 6),
            "BfaI" => array("BfaI,FspBI,MaeI,XspI", "C'TA_G", "(CTAG)", 4, 1, 2, 4),
            "BfmI" => array("BfmI,BstSFI,SfcI", "C'TRYA_G", "(CTACAG|CTATAG|CTGCAG|CTGTAG)", 6, 1, 4, 6),
            "BfrBI" => array("BfrBI", "ATG'CAT", "(ATGCAT)", 6, 3, 0, 6),
            "BfuCI" => array("BfuCI,Bsp143I,BssMI,BstMBI,DpnII,Kzo9I,MboI,NdeII,Sau3AI", "'GATC_", "(GATC)", 4, 0, 4, 4),
            "BglI" => array("BglI", "GCCN_NNN'NGGC", "(GCC.....GGC)", 11, 7, -3, 6),
            "BglII" => array("BglII", "A'GATC_T", "(AGATCT)", 6, 1, 4, 6),
            "BisI" => array("BisI,BlsI,Fnu4HI,Fsp4HI,GluI,ItaI,SatI", "GC'N_GC", "(GC.GC)", 5, 2, 1, 4),
            "BlpI" => array("BlpI,Bpu1102I,Bsp1720I,CelII", "GC'TNA_GC", "(GCT.AGC)", 7, 2, 3, 6),
            "Bme1390I" => array("Bme1390I,MspR9I,ScrFI", "CC'N_GG", "(CC.GG)", 5, 2, 1, 4),
            "Bme1580I" => array("Bme1580I,BseSI,BstSLI", "G_KGCM'C", "(GGGCAC|GGGCCC|GTGCAC|GTGCCC)", 6, 5, -4, 6),
            "BmiI" => array("BmiI,BspLI,NlaIV,PspN4I", "GGN'NCC", "(GG..CC)", 6, 3, 0, 4),
            "BmrFI" => array("BmrFI,BssKI,BstSCI,StyD4I", "'CCNGG_", "(CC.GG)", 5, 0, 5, 4),
            "BmtI" => array("BmtI", "G_CTAG'C", "(GCTAGC)", 6, 5, -4, 6),
            "BoxI" => array("BoxI,PshAI,BstPAI", "GACNN'NNGTC", "(GAC....GTC)", 10, 5, 0, 6),
            "Bpu10I" => array("Bpu10I", "CC'TNA_GC", "(CCT.AGC|GCT.AGG)", 7, 2, 3, 6),
            "Bpu14I" => array("Bpu14I,Bsp119I,BspT104I,BstBI,Csp45I,NspV,SfuI", "TT'CG_AA", "(TTCGAA)", 6, 2, 2, 6),
            "BpvUI" => array("BpvUI,MvrI,PvuI,Ple19I", "CG_AT'CG", "(CGATCG)", 6, 4, -2, 6),
            "BsaAI" => array("BsaAI,BstBAI,Ppu21I", "YAC'GTR", "(CACGTA|CACGTG|TACGTA|TACGTG)", 6, 3, 0, 6),
            "BsaBI" => array("BsaBI,Bse8I,BseJI,MamI", "GATNN'NNATC", "(GAT....ATC)", 10, 5, 0, 6),
            "BsaJI" => array("BsaJI,BseDI,BssECI", "C'CNNG_G", "(CC..GG)", 6, 1, 4, 4),
            "BsaWI" => array("BsaWI", "W'CCGG_W", "(ACCGGA|ACCGGT|TCCGGA|TCCGGT)", 6, 1, 4, 6),
            "Bse118I" => array("Bse118I,BsrFI,BssAI,Cfr10I", "R'CCGG_Y", "(ACCGGC|ACCGGT|GCCGGC|GCCGGT)", 6, 1, 4, 6),
            "BseBI" => array("BseBI,BstNI,BstOI,Bst2UI,MvaI", "CC'W_GG", "(CCAGG|CCTGG)", 5, 2, 1, 5),
            "BsePI" => array("BsePI,BssHII,PauI", "G'CGCG_C", "(GCGCGC)", 6, 1, 4, 6),
            "BseX3I" => array("BseX3I,BstZI,EagI,EclXI,Eco52I", "C'GGCC_G", "(CGGCCG)", 6, 1, 4, 6),
            "BseYI" => array("BseYI", "C'CCAG_C", "(CCCAGC|GCTGGG)", 6, 1, 4, 6),
            "Bsh1285I" => array("Bsh1285I,BsiEI,BstMCI", "CG_RY'CG", "(CGACCG|CGATCG|CGGCCG|CGGTCG)", 6, 4, -2, 6),
            "BshFI" => array("BshFI,BsnI,BspANI,BsuRI,HaeIII,PhoI", "GG'CC", "(GGCC)", 4, 2, 0, 4),
            "BsiSI" => array("BsiSI,HapII,HpaII,MspI", "C'CG_G", "(CCGG)", 4, 1, 2, 4),
            "BsiWI" => array("BsiWI,Pfl23II,PspLI", "C'GTAC_G", "(CGTACG)", 6, 1, 4, 6),
            "Bsp120I" => array("Bsp120I,PspOMI", "G'GGCC_C", "(GGGCCC)", 6, 1, 4, 6),
            "Bsp1286I" => array("Bsp1286I,MhlI,SduI", "G_DGCH'C", "(GAGCAC|GAGCTC|GAGCCC|GTGCAC|GTGCTC|GTGCCC|GGGCAC|GGGCTC|GGGCCC)", 6, 5, -4, 6),
            "Bsp1407I" => array("Bsp1407I,BsrGI,BstAUI", "T'GTAC_A", "(TGTACA)", 6, 1, 4, 6),
            "Bsp19I" => array("Bsp19I,NcoI", "C'CATG_G", "(CCATGG)", 6, 1, 4, 6),
            "Bsp68I" => array("Bsp68I,BtuMI,NruI", "TCG'CGA", "(TCGCGA)", 6, 3, 0, 6),
            "BspHI" => array("BspHI,PagI,RcaI", "T'CATG_A", "(TCATGA)", 6, 1, 4, 6),
            "BspLU11I" => array("BspLU11I,PciI,PscI", "A'CATG_T", "(ACATGT)", 6, 1, 4, 6),
            "BspMAI" => array("BspMAI,PstI", "C_TGCA'G", "(CTGCAG)", 6, 5, -4, 6),
            "BssNAI" => array("BssNAI,Bst1107I,BstZ17I", "GTA'TAC", "(GTATAC)", 6, 3, 0, 6),
            "BssSI" => array("BssSI,Bst2BI", "C'ACGA_G or C'TCGT_G", "(CACGAG|CTCGTG)", 6, 1, 4, 6),
            "BssT1I" => array("BssT1I,StyI,Eco130I,EcoT14I,ErhI", "C'CWWG_G", "(CCAAGG|CCATGG|CCTAGG|CCTTGG)", 6, 1, 4, 6),
            "Bst4CI" => array("Bst4CI,HpyCH4III,TaaI", "AC_N'GT", "(AC.GT)", 5, 3, -1, 4),
            "BstH2I" => array("BstH2I,HaeII", "R_GCGC'Y", "(AGCGCC|AGCGCT|GGCGCC|GGCGCT)", 6, 5, -4, 6),
            "BstAPI" => array("BstAPI", "GCAN_NNN'NTGC", "(GCA.....TGC)", 11, 7, -3, 6),
            "BstC8I" => array("BstC8I,Cac8I", "GCN'NGC", "(GC..GC)", 6, 3, 0, 4),
            "BstDEI" => array("BstDEI,DdeI,HpyF3I", "C'TNA_G", "(CT.AG)", 5, 1, 3, 4),
            "BstDSI" => array("BstDSI,BtgI", "C'CRYG_G", "(CCACGG|CCATGG|CCGCGG|CCGTGG)", 6, 1, 4, 6),
            "BstEII" => array("BstEII,BstPI,Eco91I,EcoO65I,PspEI", "G'GTNAC_C", "(GGT.ACC)", 7, 1, 5, 6),
            "BstENI" => array("BstENI,EcoNI,XagI", "CCTNN'N_NNAGG", "(CCT.....AGG)", 11, 5, 1, 6),
            "BstKTI" => array("BstKTI", "G_AT'C", "(GATC)", 4, 3, 2, 4),
            "BstMWI" => array("BstMWI,MwoI", "GCNN_NNN'NNGC", "(GC.......GC)", 11, 7, -3, 4),
            "BstNSI" => array("BstNSI,NspI,XceI", "R_CATG'Y", "(ACATGC|ACATGT|GCATGC|GCATGT)", 6, 5, -4, 6),
            "BstSNI" => array("BstSNI,Eco105I,SnaBI", "TAC'GTA", "(TACGTA)", 6, 3, 0, 6),
            "BstX2I" => array("BstX2I,BstYI,MflI,PsuI,XhoII", "R'GATC_Y", "(AGATCC|AGATCT|GGATCC|GGATCT)", 6, 1, 4, 6),
            "BstXI" => array("BstXI", "CCAN_NNNN'NTGG", "(CCA......TGG)", 12, 8, -4, 6),
            "CciNI" => array("CciNI,NotI", "GC'GGCC_GC", "(GCGGCCGC)", 8, 2, 4, 8),
            "Cfr42I" => array("Cfr42I,KspI,SacII,Sfr303I,SgrBI,SstII", "CC_GC'GG", "(CCGCGG)", 6, 4, -2, 6),
            "Cfr9I" => array("Cfr9I,TspMI,XmaI,XmaCI", "C'CCGG_G", "(CCCGGG)", 6, 1, 4, 6),
            "CpoI" => array("CpoI,CspI,RsrII,Rsr2I", "CG'GWC_CG", "(CGGACCG|CGGTCCG)", 7, 2, 3, 7),
            "Csp6I" => array("Csp6I,CviQI,RsaNI", "G'TA_C", "(GTAC)", 4, 1, 2, 4),
            "CviAII" => array("CviAII,FaeI,Hin1II,Hsp92II,NlaIII", "_CATG'", "(CATG)", 4, 4, -4, 4),
            "CviJI" => array("CviJI,CviKI-1", "RG'CY", "(AGCC|AGCT|GGCC|GGCT)", 4, 2, 0, 4),
            "DinI" => array("DinI,Mly113I,NarI", "GG'CG_CC", "(GGCGCC)", 6, 2, 2, 6),
            "DpnI" => array("DpnI,MalI", "GA'TC", "(GATC)", 4, 2, 0, 4),
            "DraI" => array("DraI", "TTT'AAA", "(TTTAAA)", 6, 3, 0, 6),
            "Ecl136II" => array("Ecl136II,EcoICRI", "GAG'CTC", "(GAGCTC)", 6, 3, 0, 6),
            "Eco32I" => array("Eco32I,EcoRV", "GAT'ATC", "(GATATC)", 6, 3, 0, 6),
            "EcoO109I" => array("EcoO109I,DraII", "RG'GNC_CY", "(AGG.CCC|AGG.CCT|GGG.CCC|GGG.CCT)", 7, 2, 3, 6),
            "EcoRI" => array("EcoRI", "G'AATT_C", "(GAATTC)", 6, 1, 4, 6),
            "EcoT22I" => array("EcoT22I,Mph1103I,NsiI,Zsp2I", "A_TGCA'T", "(ATGCAT)", 6, 5, -4, 6),
            "EgeI" => array("EgeI,EheI,SfoI", "GGC'GCC", "(GGCGCC)", 6, 3, 0, 6),
            "FatI" => array("FatI", "'CATG_", "(CATG)", 4, 0, 4, 4),
            "FauNDI" => array("FauNDI,NdeI", "CA'TA_TG", "(CATATG)", 6, 2, 2, 6),
            "FseI" => array("FseI,RigI", "GG_CCGG'CC", "(GGCCGGCC)", 8, 6, -4, 8),
            "FspAI" => array("FspAI", "RTGC'GCAY", "(ATGCGCAC|ATGCGCAT|GTGCGCAC|GTGCGCAT)", 8, 4, 0, 8),
            "GlaI" => array("GlaI", "GC'GC", "(GCGC)", 4, 2, 0, 4),
            "Hin6I" => array("Hin6I,HinP1I,HspAI", "G'CG_C", "(GCGC)", 4, 1, 2, 4),
            "HincII" => array("HincII,HindII", "GTY'RAC", "(GTCAAC|GTCGAC|GTTAAC|GTTGAC)", 6, 3, 0, 6),
            "HindIII" => array("HindIII", "A'AGCT_T", "(AAGCTT)", 6, 1, 4, 6),
            "HinfI" => array("HinfI", "G'ANT_C", "(GA.TC)", 5, 1, 3, 4),
            "HpaI" => array("HpaI,KspAI", "GTT'AAC", "(GTTAAC)", 6, 3, 0, 6),
            "Hpy188I" => array("Hpy188I", "TC_N'GA", "(TC.GA)", 5, 3, -1, 4),
            "Hpy188III" => array("Hpy188III", "TC'NN_GA", "(TC..GA)", 6, 2, 2, 4),
            "Hpy8I" => array("Hpy8I", "GTN'NAC", "(GT..AC)", 6, 3, 0, 4),
            "Hpy99I" => array("Hpy99I", "_CGWCG'", "(CGACG|CGTCG)", 5, 5, -5, 5),
            "HpyCH4IV" => array("HpyCH4IV,MaeII", "A'CG_T", "(ACGT)", 4, 1, 2, 4),
            "HpyCH4V" => array("HpyCH4V", "TG'CA", "(TGCA)", 4, 2, 0, 4),
            "HpyF10VI" => array("HpyF10VI", "GCNN_NNN'NNGC", "(GC.......GC)", 11, 7, -3, 4),
            "KasI" => array("KasI", "G'GCGC_C", "(GGCGCC)", 6, 1, 4, 6),
            "KpnI" => array("KpnI", "G_GTAC'C", "(GGTACC)", 6, 5, -4, 6),
            "MabI" => array("MabI,SexAI", "A'CCWGG_T", "(ACCAGGT|ACCTGGT)", 7, 1, 5, 7),
            "MaeIII" => array("MaeIII", "'GTNAC_", "(GT.AC)", 5, 0, 5, 4),
            "MfeI" => array("MfeI,MunI", "C'AATT_G", "(CAATTG)", 6, 1, 4, 6),
            "MluI" => array("MluI", "A'CGCG_T", "(ACGCGT)", 6, 1, 4, 6),
            "MreI" => array("MreI", "CG'CCGG_CG", "(CGCCGGCG)", 8, 2, 4, 8),
            "MroNI" => array("MroNI,NgoMIV", "G'CCGG_C", "(GCCGGC)", 6, 1, 4, 6),
            "MseI" => array("MseI,Tru1I,Tru9I", "T'TA_A", "(TTAA)", 4, 1, 2, 4),
            "MslI" => array("MslI,RseI,SmiMI", "CAYNN'NNRTG", "(CAC....ATG|CAC....GTG|CAT....ATG|CAT....GTG)", 10, 5, 0, 6),
            "MspA1I" => array("MspA1I", "CMG'CKG", "(CAGCGG|CAGCTG|CCGCGG|CCGCTG)", 6, 3, 0, 6),
            "MssI" => array("MssI,PmeI", "GTTT'AAAC", "(GTTTAAAC)", 8, 4, 0, 8),
            "NaeI" => array("NaeI,PdiI", "GCC'GGC", "(GCCGGC)", 6, 3, 0, 6),
            "NmuCI" => array("NmuCI,Tsp45I", "'GTSAC_", "(GTCAC|GTGAC)", 5, 0, 5, 5),
            "PacI" => array("PacI", "TTA_AT'TAA", "(TTAATTAA)", 8, 5, -2, 8),
            "PaeR7I" => array("PaeR7I,Sfr274I,SlaI,StrI,TliI,XhoI", "C'TCGA_G", "(CTCGAG)", 6, 1, 4, 6),
            "PasI" => array("PasI", "CC'CWG_GG", "(CCCAGGG|CCCTGGG)", 7, 2, 3, 7),
            "PfeI" => array("TfiI,PfeI", "G'AWT_C", "(GAATC|GATTC)", 5, 1, 3, 5),
            "PfoI" => array("PfoI", "T'CCNGG_A", "(TCC.GGA)", 7, 1, 5, 6),
            "PpuMI" => array("PpuMI,Psp5II,PspPPI", "RG'GWC_CY", "(AGGACCC|AGGACCT|AGGTCCC|AGGTCCT|GGGACCC|GGGACCT|GGGTCCC|GGGTCCT)", 7, 2, 3, 7),
            "PsiI" => array("PsiI", "TTA'TAA", "(TTATAA)", 6, 3, 0, 6),
            "Psp124BI" => array("Psp124BI,SacI,SstI", "G_AGCT'C", "(GAGCTC)", 6, 5, -44, 6),
            "PspXI" => array("PspXI", "VC'TCGA_GB", "(ACTCGAGC|ACTCGAGG|ACTCGAGT|CCTCGAGC|CCTCGAGG|CCTCGAGT|GCTCGAGC|GCTCGAGG|GCTCGAGT)", 8, 2, 4, 8),
            "PvuII" => array("PvuII", "CAG'CTG", "(CAGCTG)", 6, 3, 0, 6),
            "SalI" => array("SalI", "G'TCGA_C", "(GTCGAC)", 6, 1, 4, 6),
            "SanDI" => array("SanDI", "GG'GWC_CC", "(GGGACCC|GGGTCCC)", 7, 2, 3, 7),
            "SbfI" => array("SbfI,SdaI,Sse8387I", "CC_TGCA'GG", "(CCTGCAGG)", 8, 6, -4, 8),
            "SetI" => array("SetI", "_ASST'", "(AGGT|AGCT|ACGT|ACCT)", 4, 4, -4, 4),
            "SfiI" => array("SfiI", "GGCCN_NNN'NGGCC", "(GGCC.....GGCC)", 13, 8, -3, 8),
            "SgrAI" => array("SgrAI", "CR'CCGG_YG", "(CACCGGCG|CACCGGTG|CGCCGGCG|CGCCGGTG)", 8, 2, 4, 8),
            "SgrDI" => array("SgrDI", "CG'TCGA_CG", "(CGTCGACG)", 8, 2, 4, 8),
            "SmaI" => array("SmaI", "CCC'GGG", "(CCCGGG)", 6, 3, 0, 6),
            "SmiI" => array("SmiI,SwaI", "ATTT'AAAT", "(ATTTAAAT)", 8, 4, 0, 8),
            "SmlI" => array("SmlI,SmoI", "C'TYRA_G", "(CTCAAG|CTCGAG|CTTAAG|CTTGAG)", 6, 1, 4, 6),
            "SrfI" => array("SrfI", "GCCC'GGGC", "(GCCCGGGC)", 8, 4, 0, 8),
            "Sse9I" => array("Sse9I,TasI,Tsp509I,TspEI", "'AATT_", "(AATT)", 4, 0, 4, 4),
            "SspI" => array("SspI", "AAT'ATT", "(AATATT)", 6, 3, 0, 6),
            "TaiI" => array("TaiI", "_ACGT'", "(ACGT)", 4, 4, -4, 4),
            "TaqI" => array("TaqI", "T'CG_A", "(TCGA)", 4, 1, 2, 4),
            "TatI" => array("TatI", "W'GTAC_W", "(AGTACA|AGTACT|TGTACA|TGTACT)", 6, 1, 4, 6),
            "TauI" => array("TauI", "G_CSG'C", "(GCCGC|GCGGC)", 5, 4, -3, 5),
            "TspRI" => array("TspRI", "_NNCASTGNN'", "(..CACTG..|..CAGTG..)", 9, 9, -9, 5),
            "XbaI" => array("XbaI", "T'CTAG_A", "(TCTAGA)", 6, 1, 4, 6),
            "XcmI" => array("XcmI", "CCANNNN_N'NNNNTGG", "(CCA.........TGG)", 15, 8, -1, 6),
            "ZraI" => array("ZraI", "GAC'GTC", "(GACGTC)", 6, 3, 0, 6),
        );
        return $enzymes_array;
    }

    private function getEnzymesTypeIIs() {
        // Two lines for each endonuclease
        $enzymes_array = Array(
            "AarI" => array("AarI", "CACCTGCNNNN'NNNN_", "(CACCTGC........)", 15, 11, 4, 7),
            "AarI@" => array("", "", "(........GCAGGTG)", 15, 0, 4, 7),
            "Acc36I" => array("Acc36I,BfuAI,BspMI,BveI", "ACCTGCNNNN'NNNN_", "(ACCTGC........)", 14, 10, 4, 6),
            "Acc36I@" => array("", "", "(........GCAGGT)", 14, 0, 4, 6),
            "AclWI" => array("AclWI,AlwI,BspPI", "GGATCNNNN'N_", "(GGATC.....)", 10, 9, 1, 5),
            "AclWI@" => array("", "", "(.....GATCC)", 10, 0, 1, 5),
            "AcuI" => array("AcuI,Eco57I", "R'AATT_Y", "(CTGAAGNNNNNNNNNNNNNN_NN')", 22, 22, -2, 6),
            "AcuI@" => array("", "", "(................CTTCAG)", 22, 2, -2, 6),
            "Alw26I" => array("Alw26I,BsmAI,BstMAI", "GTCTCN'NNNN_", "(GTCTC.....)", 10, 6, 4, 5),
            "Alw26I@" => array("", "", "(.....GAGAC)", 10, 0, 4, 5),
            "AsuHPI" => array("AsuHPI,HphI", "GTGANNNNNNN_N'", "(GGTGA........)", 13, 13, -1, 5),
            "AsuHPI@" => array("", "", "(........TCACC)", 13, 1, -1, 5),
            "BbsI" => array("BbsI,BpiI,BpuAI,BstV2I", "GAAGACNN'NNNN_", "(GAAGAC......)", 12, 8, 4, 6),
            "BbsI@" => array("", "", "(......GTCTTC)", 12, 0, 4, 6),
            "BbvI" => array("BbvI,BseXI,BstV1I", "GCAGCNNNNNNNN'NNNN_", "(GCAGC............)", 17, 13, 4, 5),
            "BbvI@" => array("", "", "(............GCTGC)", 17, 0, 4, 5),
            "BccI" => array("BccI", "CCATCNNNN'N_", "(CCATC.....)", 10, 9, 1, 5),
            "BccI@" => array("", "", "(.....GATGG)", 10, 0, 1, 5),
            "BceAI" => array("BceAI", "ACGGCNNNNNNNNNNNN'NN_", "(ACGGC..............)", 19, 17, 2, 5),
            "BceAI@" => array("", "", "(..............GCCGT)", 19, 0, 2, 5),
            "BciVI" => array("BciVI,BfuI", "GTATCCNNNNN_N'", "(GTATCC......)", 12, 12, -1, 6),
            "BciVI@" => array("", "", "(......GGATAC)", 12, 1, -1, 6),
            "BfiI" => array("BfiI,BmrI,BmuI", "ACTGGGNNNN_N'", "(ACTGGG.....)", 11, 11, -1, 6),
            "BfiI@" => array("", "", "(.....CCCAGT)", 11, 1, -1, 6),
            "BpmI" => array("BpmI,GsuI", "CTGGAGNNNNNNNNNNNNNN_NN'", "(CTGGAG................)", 22, 22, -2, 6),
            "BpmI@" => array("", "", "(................CTCCAG)", 22, 2, -2, 6),
            "BpuEI" => array("BpuEI", "CTTGAGNNNNNNNNNNNNNN_NN'", "(CTTGAG................)", 22, 22, -2, 6),
            "BpuEI@" => array("", "", "(................CTCAAG)", 22, 2, -2, 6),
            "BsaI" => array("BsaI,Bso31I,BspTNI,Eco31I", "GGTCTCN'NNNN_", "(GGTCTC.....)", 11, 7, 4, 6),
            "BsaI@" => array("", "", "(.....GAGACC)", 11, 0, 4, 6),
            "BsaMI" => array("BsaMI,BsmI,Mva1269I,PctI", "GAATG_CN'", "(GAATGC.)", 7, 7, -2, 6),
            "BsaMI@" => array("", "", "(.GCATTC)", 7, 2, -2, 6),
            "Bse1I" => array("Bse1I,BseNI,BsrI,BsrSI", "ACTG_GN'", "(ACTGG.)", 6, 6, -2, 5),
            "Bse1I@" => array("", "", "(.CCAGT)", 6, 2, -2, 5),
            "Bse3DI" => array("Bse3DI,BseMI,BsrDI", "GCAATG_NN'", "(GCAATG..)", 8, 8, -2, 6),
            "Bse3DI@" => array("", "", "(..CATTGC)", 8, 2, -2, 6),
            "BseGI" => array("BseGI,BstF5I", "GGATG_NN'", "(GGATG..)", 7, 7, -2, 5),
            "BseGI@" => array("", "", "(..CATCC)", 7, 2, -2, 5),
            "BseMII" => array("BseMII", "CTCAGNNNNNNNN_NN'", "(CTCAG..........)", 15, 15, -2, 5),
            "BseMII@" => array("", "", "(..........CTGAG)", 15, 2, -2, 5),
            "BseRI" => array("BseRI", "GAGGAGNNNNNNNN_NN'", "(GAGGAG..........)", 16, 16, -2, 6),
            "BseRI@" => array("", "", "(..........CTCCTC)", 16, 2, -2, 6),
            "BsgI" => array("BsgI", "GTGCAGNNNNNNNNNNNNNN_NN'", "(GTGCAG................)", 22, 22, -2, 6),
            "BsgI@" => array("", "", "(................CTGCAC)", 22, 2, -2, 6),
            "BslFI" => array("BslFI,FaqI", "GGGACNNNNNNNNNN'NNNN_", "(GGGAC..............)", 19, 15, 4, 5),
            "BslFI@" => array("", "", "(..............GTCCC)", 19, 0, 4, 5),
            "BsmBI" => array("BsmBI,Esp3I", "CGTCTCN'NNNN_", "(CGTCTC.....)", 11, 7, 4, 6),
            "BsmBI@" => array("", "", "(.....GAGACG)", 11, 0, 4, 6),
            "BsmFI" => array("BsmFI", "GGGACNNNNNNNNNN'NNNN_", "(GGGAC..............)", 19, 15, 4, 5),
            "BsmFI@" => array("", "", "(..............GTCCC)", 19, 0, 4, 5),
            "BspCNI" => array("BspCNI", "CTCAGNNNNNNN_NN'", "(CTCAG.........)", 14, 14, -2, 5),
            "BspCNI@" => array("", "", "(.........CTGAG)", 14, 2, -2, 5),
            "BspQI" => array("BspQI,LguI,PciSI,SapI", "GCTCTTCN'NNN_", "(GCTCTTC....)", 11, 8, 3, 7),
            "BspQI@" => array("", "", "(....GAAGAGC)", 11, 0, 3, 7),
            "Bst6I" => array("Bst6I,Eam1104I,EarI,Ksp632I", "CTCTTCN'NNN_", "(CTCTTC....)", 10, 7, 3, 6),
            "Bst6I@" => array("", "", "(....GAAGAG)", 10, 0, 3, 6),
            "BtgZI" => array("BtgZI", "GCGATGNNNNNNNNNN'NNNN_", "(GCGATG..............)", 20, 16, 4, 6),
            "BtgZI@" => array("", "", "(..............CATCGC)", 20, 0, 4, 6),
            "BtsI" => array("BtsI", "GCAGTG_NN'", "(GCAGTG..)", 8, 8, -2, 6),
            "BtsI@" => array("", "", "(..CACTGC)", 8, 2, -2, 6),
            "EciI" => array("EciI", "GGCGGANNNNNNNNN_NN'", "(GGCGGA...........)", 17, 17, -2, 6),
            "EciI@" => array("", "", "(...........TCCGCC)", 17, 2, -2, 6),
            "Eco57MI" => array("Eco57MI", "CTGRAGNNNNNNNNNNNNNN_NN'", "(CTGAAG................|CTGGAG................)", 22, 22, -2, 6),
            "Eco57MI@" => array("", "", "(................CTCCAG|................CTTCAG)", 22, 2, -2, 6),
            "EcoP15I" => array("EcoP15I", "CAGCAGNNNNNNNNNNNNNNNNNNNNNNNNN'NN_", "(CAGCAG...........................)", 33, 31, 2, 6),
            "EcoP15I@" => array("", "", "(...........................CTGCTG)", 33, 0, 2, 6),
            "FauI" => array("FauI,SmuI", "CCCGCNNNN'NN_", "(CCCGC......)", 11, 9, 2, 5),
            "FauI@" => array("", "", "(......GCGGG)", 11, 0, 2, 5),
            "FokI" => array("FokI,BtsCI", "GGATGNNNNNNNNN'NNNN_", "(GGATG.............)", 18, 14, 4, 5),
            "FokI@" => array("", "", "(.............CATCC)", 18, 0, 4, 5),
            "HgaI" => array("HgaI,CseI", "GACGCNNNNN'NNNNN_", "(GACGC..........)", 15, 10, 5, 5),
            "HgaI@" => array("", "", "(..........GCGTC)", 15, 0, 5, 5),
            "HpyAV" => array("HpyAV", "CCTTCNNNNN_N'", "(GACGC..........)", 11, 11, -1, 5),
            "HpyAV@" => array("", "", "(......GAAGG)", 11, 1, -1, 5),
            "LweI" => array("LweI,SfaNI", "GCATCNNNNN'NNNN_", "(GCATC.........)", 14, 10, 4, 5),
            "LweI@" => array("", "", "(.........GATGC)", 14, 0, 4, 5),
            "MboII" => array("MboII", "GAAGANNNNNNN_N'", "(GAAGA........)", 13, 13, -1, 5),
            "MboII@" => array("", "", "(........TCTTC)", 13, 1, -1, 5),
            "MlyI" => array("MlyI,SchI", "GAGTCNNNNN'", "(GAGTC.....)", 10, 10, 0, 5),
            "MlyI@" => array("", "", "(.....GACTC)", 10, 0, 0, 5),
            "MmeAIII" => array("MmeAIII", "TCCRACNNNNNNNNNNNNNNNNNN_NN'", "(GCCGAG.....................)", 27, 27, -2, 6),
            "MmeAIII@" => array("", "", "(.....................CTCGGC)", 27, 2, -2, 6),
            "MmeI" => array("MmeI", "TCCRACNNNNNNNNNNNNNNNNNN_NN'", "(TCCAAC....................|TCCGAC....................)", 26, 26, -2, 6),
            "MmeI@" => array("", "", "(....................GTCGGA|....................GTTGGA)", 26, 2, -2, 6),
            "MnlI" => array("MnlI", "CCTCNNNNNN_N'", "(CCTC.......)", 11, 11, -1, 4),
            "MnlI@" => array("", "", "(.......GAGG)", 11, 1, -1, 4),
            "PleI" => array("PleI,PpsI", "GAGTCNNNN'N_", "(GAGTC.....)", 10, 9, 1, 5),
            "PleI@" => array("", "", "(.....GACTC)", 10, 0, 1, 5),
            "TaqII" => array("TaqII", "GACCGANNNNNNNNN_NN' or CACCCANNNNNNNNN_NN'", "(GACCGA...........|CACCCA...........)", 17, 17, -2, 6),
            "TaqII@" => array("", "", "(...........TCGGTC|...........TGGGTG)", 17, 2, -2, 6),
            "TsoI" => array("TsoI", "TARCCANNNNNNNNN_NN'", "(TAACCA...........|TAGCCA...........)", 17, 17, -2, 6),
            "TsoI@" => array("", "", "(...........TGGTTA|...........TGGCTA)", 17, 2, -2, 6),
            "TspDTI" => array("TspDTI", "ATGAANNNNNNNNN_NN'", "(ATGAA...........)", 16, 16, -2, 5),
            "TspDTI@" => array("", "", "(...........TTCAT)", 16, 2, -2, 5),
            "TspGWI" => array("TspGWI", "ACGGANNNNNNNNN_NN'", "(ACGGA...........)", 16, 16, -2, 5),
            "TspGWI@" => array("", "", "(...........TCCGT)", 16, 2, -2, 5),
        );
        return $enzymes_array;
    }

    private function getEnzymesTypeIIb() {
        $enzymes_array = Array(
            "AjuI#" => array("AjuI", "_NNNNN'NNNNNNNGAANNNNNNNTTGGNNNNNN_NNNNN_'\n   Generates two cuts", "(............GAA.......TTGG...........|...........CCAA.......TTC............)", 37, 5, -5, 7),
            "AlfI#" => array("AlfI", "_NN'NNNNNNNNNNCGANNNNNNTGCNNNNNNNNNN_NN'\n   Generates two cuts", "(............CGA......TGC............|............GCA......TCG............)", 36, 2, -2, 6),
            "AloI#" => array("AloI", "_NNNNN'NNNNNNNGAACNNNNNNTCCNNNNNNN_NNNNN'\n   Generates two cuts", "(............GAAC......TCC............|............GGA......GTTC............)", 33, 5, -5, 7),
            "BarI#" => array("BarI", "_NNNNN'NNNNNNNGAAGNNNNNNTACNNNNNNN_NNNNN'\n   Generates two cuts", "(............GAAG......TAC............|............GTA......CTTC............)", 37, 5, -5, 7),
            "BaeI#" => array("BaeI", "_NNNNN'NNNNNNNNNNACNNNNGTAYCNNNNNNN_NNNNN'\n   Generates two cuts", "(...............AC....GTACC............|...............AC....GTATC............|............GATAC....GT...............|............GGTAC....GT...............)", 38, 5, -5, 7),
            "BcgI#" => array("BcgI", "_NN'NNNNNNNNNNCGANNNNNNTGCNNNNNNNNNN_NN'\n   Generates two cuts", "(............CGA......TGC............|............GCA......TCG............)", 36, 2, -2, 6),
            "BdaI#" => array("BdaI", "_NN'NNNNNNNNNNTGANNNNNNTCANNNNNNNNNN_NN'\n   Generates two cuts", "(............TGA......TCA............)", 36, 2, -2, 6),
            "BsaXI#" => array("BsaXI", "_NNN'NNNNNNNNNACNNNNNCTCCNNNNNNN_NNN'\n   Generates two cuts", "(............AC.....CTCC..........|..........GGAG.....GT............)", 33, 3, -3, 6),
            "CspCI#" => array("CspCI", "_NN'NNNNNNNNNNNCAANNNNNGTGGNNNNNNNNNN_NN'\n   Generates two cuts", "(.............CAA.....GTGG............|............GCA.....TCG.............)", 37, 2, -2, 7),
            "Hin4I#" => array("Hin4I", "_NNNNN'NNNNNNNNGAYNNNNNVTCNNNNNNNN_NNNNN'\n   Generates two cuts", "(.............GAC.....ATC.............|.............GAC.....CTC.............|.............GAC.....GTC.............|.............GAT.....ATC.............|.............GAT.....CTC.............|.............GAT.....GTC.............|.............GAG.....ATC.............|.............GAG.....ATC.............)", 37, 5, -5, 6),
            "PpiI#" => array("PpiI", "_NNNNN'NNNNNNNGAACNNNNNCTCNNNNNNNN_NNNNN'\n   Generates two cuts", "(............GAAC.....CTC.............|.............GAG.....GTTC............)", 37, 5, -5, 7),
            "PsrI#" => array("PsrI", "_NNNNN'NNNNNNNGAACNNNNNNTACNNNNNNN_NNNNN'\n   Generates two cuts", "(............GAAC......TAC............|............GTA......GTTC............)", 35, 5, -5, 7),
            "TstI#" => array("TstI", "_NNNNN'NNNNNNNNCACNNNNNNTCCNNNNNNN_NNNNN'\n   Generates two cuts", "(.............CAC......TCC............|............GGA......GTG.............)", 37, 5, -5, 7),
        );

        return $enzymes_array;
    }

    private function endonuclease_vendors() {
        // This function return the list of sellers for all endonucleases included in this script

        $vendors = array(
            "AarI" => "F",
            "AasI" => "F",
            "AatI" => "O",
            "AatII" => "AFGIKMNORV",
            "AbsI" => "I",
            "AccI" => "ABGJKMNORSUWX",
            "AccII" => "AJK",
            "AccIII" => "GJKRW",
            "Acc16I" => "IV",
            "Acc36I" => "I",
            "Acc65I" => "FGINRVW",
            "AccB1I" => "IV",
            "AccB7I" => "IRV",
            "AccBSI" => "IV",
            "AciI" => "N",
            "AclI" => "INV",
            "AclWI" => "I",
            "AcoI" => "I",
            "AcsI" => "IMV",
            "AcuI" => "IN",
            "AcvI" => "QX",
            "AcyI" => "JM",
            "AdeI" => "F",
            "AfaI" => "AK",
            "AfeI" => "IN",
            "AfiI" => "V",
            "AflII" => "AJKNO",
            "AflIII" => "GMNSW",
            "AgeI" => "JNR",
            "AhdI" => "N",
            "AhlI" => "IV",
            "AjiI" => "F",
            "AjnI" => "I",
            "AjuI" => "F",
            "AleI" => "N",
            "AlfI" => "F",
            "AloI" => "F",
            "AluI" => "ABFGHIJKMNOQRSUVWXY",
            "AluBI" => "I",
            "AlwI" => "N",
            "Alw21I" => "F",
            "Alw26I" => "FR",
            "Alw44I" => "FJMORS",
            "AlwNI" => "N",
            "Ama87I" => "IV",
            "Aor13HI" => "K",
            "Aor51HI" => "AK",
            "ApaI" => "ABFGIJKMNOQRSUVWX",
            "ApaLI" => "AKNU",
            "ApeKI" => "N",
            "ApoI" => "N",
            "AscI" => "GNW",
            "AseI" => "JNO",
            "AsiGI" => "IV",
            "AsiSI" => "N",
            "AspI" => "M",
            "Asp700I" => "M",
            "Asp718I" => "M",
            "AspA2I" => "IV",
            "AspEI" => "M",
            "AspLEI" => "IV",
            "AspS9I" => "IV",
            "AssI" => "U",
            "AsuC2I" => "I",
            "AsuHPI" => "IV",
            "AsuNHI" => "IV",
            "AvaI" => "ABGJMNORSUWX",
            "AvaII" => "AGJKMNRSWY",
            "AviII" => "M",
            "AvrII" => "N",
            "AxyI" => "J",
            "BaeI" => "N",
            "BalI" => "AJKR",
            "BamHI" => "ABFGHIJKMNOQRSUVWXY",
            "BanI" => "NORU",
            "BanII" => "AGKMNOQRSWX",
            "BanIII" => "O",
            "BarI" => "I",
            "BasI" => "U",
            "BauI" => "F",
            "BbeI" => "AK",
            "BbrPI" => "MO",
            "BbsI" => "N",
            "BbuI" => "R",
            "BbvI" => "N",
            "Bbv12I" => "IV",
            "BbvCI" => "N",
            "BccI" => "N",
            "BceAI" => "N",
            "BcgI" => "N",
            "BciVI" => "N",
            "BclI" => "FGJMNORSUWY",
            "BcnI" => "FK",
            "BcuI" => "F",
            "BdaI" => "F",
            "BfaI" => "N",
            "BfiI" => "F",
            "BfmI" => "F",
            "BfrI" => "MO",
            "BfuI" => "F",
            "BfuAI" => "N",
            "BfuCI" => "N",
            "BglI" => "AFGHIJKMNOQRSUVWXY",
            "BglII" => "ABFGHIJKMNOQRSUVWXY",
            "BisI" => "I",
            "BlnI" => "AKMS",
            "BlpI" => "N",
            "BlsI" => "I",
            "BmcAI" => "V",
            "Bme18I" => "IV",
            "Bme1390I" => "F",
            "Bme1580I" => "N",
            "BmeRI" => "V",
            "BmeT110I" => "K",
            "BmgBI" => "N",
            "BmgT120I" => "K",
            "BmiI" => "V",
            "BmrI" => "N",
            "BmrFI" => "V",
            "BmtI" => "INV",
            "BmuI" => "I",
            "BoxI" => "F",
            "BpiI" => "F",
            "BplI" => "F",
            "BpmI" => "IN",
            "Bpu10I" => "FINV",
            "Bpu14I" => "IV",
            "Bpu1102I" => "AFK",
            "BpuAI" => "M",
            "BpuEI" => "N",
            "BpuMI" => "V",
            "BpvUI" => "V",
            "BsaI" => "N",
            "Bsa29I" => "I",
            "BsaAI" => "N",
            "BsaBI" => "N",
            "BsaHI" => "N",
            "BsaJI" => "N",
            "BsaMI" => "GR",
            "BsaWI" => "N",
            "BsaXI" => "N",
            "Bsc4I" => "I",
            "Bse1I" => "IV",
            "Bse8I" => "IV",
            "Bse21I" => "IV",
            "Bse118I" => "IV",
            "BseAI" => "CM",
            "BseBI" => "C",
            "BseCI" => "C",
            "BseDI" => "F",
            "Bse3DI" => "IV",
            "BseGI" => "F",
            "BseJI" => "F",
            "BseLI" => "F",
            "BseMI" => "F",
            "BseMII" => "F",
            "BseNI" => "F",
            "BsePI" => "IV",
            "BseRI" => "N",
            "BseSI" => "F",
            "BseXI" => "F",
            "BseX3I" => "IV",
            "BseYI" => "N",
            "BsgI" => "N",
            "Bsh1236I" => "F",
            "Bsh1285I" => "F",
            "BshFI" => "C",
            "BshNI" => "F",
            "BshTI" => "F",
            "BshVI" => "V",
            "BsiEI" => "N",
            "BsiHKAI" => "N",
            "BsiHKCI" => "QX",
            "BsiSI" => "C",
            "BsiWI" => "MNO",
            "BsiYI" => "M",
            "BslI" => "GNW",
            "BslFI" => "I",
            "BsmI" => "JMNOSW",
            "BsmAI" => "N",
            "BsmBI" => "N",
            "BsmFI" => "N",
            "BsnI" => "V",
            "Bso31I" => "IV",
            "BsoBI" => "N",
            "Bsp13I" => "IV",
            "Bsp19I" => "IV",
            "Bsp68I" => "F",
            "Bsp119I" => "F",
            "Bsp120I" => "F",
            "Bsp143I" => "F",
            "Bsp1286I" => "JKNR",
            "Bsp1407I" => "FK",
            "Bsp1720I" => "IV",
            "BspACI" => "I",
            "BspANI" => "X",
            "BspCNI" => "N",
            "BspDI" => "N",
            "BspEI" => "N",
            "BspFNI" => "I",
            "BspHI" => "N",
            "BspLI" => "F",
            "BspLU11I" => "M",
            "BspMI" => "N",
            "BspMAI" => "X",
            "BspOI" => "F",
            "BspPI" => "F",
            "BspQI" => "N",
            "BspTI" => "F",
            "BspT104I" => "K",
            "BspT107I" => "K",
            "BspTNI" => "QX",
            "BspXI" => "GW",
            "BsrI" => "N",
            "BsrBI" => "N",
            "BsrDI" => "N",
            "BsrFI" => "N",
            "BsrGI" => "N",
            "BsrSI" => "R",
            "BssAI" => "C",
            "BssECI" => "I",
            "BssHII" => "AJKMNOQRSX",
            "BssKI" => "N",
            "BssMI" => "V",
            "BssNI" => "V",
            "BssNAI" => "IV",
            "BssSI" => "N",
            "BssT1I" => "IV",
            "Bst6I" => "IV",
            "Bst98I" => "R",
            "Bst1107I" => "FKM",
            "BstACI" => "I",
            "BstAPI" => "IN",
            "BstAUI" => "IV",
            "BstBI" => "N",
            "Bst2BI" => "IV",
            "BstBAI" => "IV",
            "Bst4CI" => "IV",
            "BstC8I" => "I",
            "BstDEI" => "IV",
            "BstDSI" => "IV",
            "BstEII" => "GHJMNORSUW",
            "BstENI" => "IV",
            "BstF5I" => "IV",
            "BstFNI" => "IV",
            "BstH2I" => "IV",
            "BstHHI" => "IV",
            "BstKTI" => "I",
            "BstMAI" => "IV",
            "BstMBI" => "IV",
            "BstMCI" => "IV",
            "BstMWI" => "I",
            "BstNI" => "N",
            "BstNSI" => "IV",
            "BstOI" => "R",
            "BstPI" => "K",
            "BstPAI" => "IV",
            "BstSCI" => "I",
            "BstSFI" => "I",
            "BstSLI" => "I",
            "BstSNI" => "IV",
            "BstUI" => "N",
            "Bst2UI" => "IV",
            "BstV1I" => "I",
            "BstV2I" => "IV",
            "BstXI" => "AFGHIJKMNOQRVWX",
            "BstX2I" => "IV",
            "BstYI" => "N",
            "BstZI" => "R",
            "BstZ17I" => "N",
            "Bsu15I" => "F",
            "Bsu36I" => "NR",
            "BsuRI" => "FI",
            "BsuTUI" => "X",
            "BtgI" => "N",
            "BtgZI" => "N",
            "BtrI" => "IV",
            "BtsI" => "N",
            "BtsCI" => "N",
            "BtuMI" => "V",
            "BveI" => "F",
            "Cac8I" => "N",
            "CaiI" => "F",
            "CciNI" => "IV",
            "CelII" => "M",
            "CfoI" => "MRS",
            "CfrI" => "F",
            "Cfr9I" => "FO",
            "Cfr10I" => "FGKO",
            "Cfr13I" => "AFO",
            "Cfr42I" => "F",
            "ClaI" => "ABHKMNRSU",
            "CpoI" => "AFK",
            "CseI" => "F",
            "CspI" => "OR",
            "Csp6I" => "F",
            "Csp45I" => "OR",
            "CspAI" => "C",
            "CspCI" => "N",
            "CviAII" => "N",
            "CviJI" => "QX",
            "CviKI-1" => "N",
            "CviQI" => "N",
            "DdeI" => "BGMNORSW",
            "DinI" => "V",
            "DpnI" => "BEFGMNRSW",
            "DpnII" => "N",
            "DraI" => "ABFGIJKMNOQRSUVWXY",
            "DraII" => "GMW",
            "DraIII" => "GIMNVW",
            "DrdI" => "N",
            "DriI" => "I",
            "DseDI" => "IV",
            "EaeI" => "AKMN",
            "EagI" => "GNW",
            "Eam1104I" => "F",
            "Eam1105I" => "FK",
            "EarI" => "N",
            "EciI" => "N",
            "Ecl136II" => "F",
            "EclHKI" => "R",
            "EclXI" => "MS",
            "Eco24I" => "F",
            "Eco31I" => "F",
            "Eco32I" => "F",
            "Eco47I" => "FO",
            "Eco47III" => "FGMORW",
            "Eco52I" => "FKO",
            "Eco57I" => "F",
            "Eco72I" => "F",
            "Eco81I" => "AFKO",
            "Eco88I" => "F",
            "Eco91I" => "F",
            "Eco105I" => "FO",
            "Eco130I" => "F",
            "Eco147I" => "F",
            "EcoICRI" => "IRV",
            "Eco57MI" => "F",
            "EcoNI" => "N",
            "EcoO65I" => "K",
            "EcoO109I" => "AFJKN",
            "EcoP15I" => "N",
            "EcoRI" => "ABCFGHIJKMNOQRSUVWXY",
            "EcoRII" => "FJMOS",
            "EcoRV" => "ABCGHIJKMNOQRSUVWXY",
            "EcoT14I" => "K",
            "EcoT22I" => "AKO",
            "EcoT38I" => "J",
            "EgeI" => "I",
            "EheI" => "FO",
            "ErhI" => "IV",
            "Esp3I" => "F",
            "FaeI" => "I",
            "FalI" => "I",
            "FaqI" => "F",
            "FatI" => "IN",
            "FauI" => "IN",
            "FauNDI" => "IV",
            "FbaI" => "AK",
            "FblI" => "IV",
            "Fnu4HI" => "N",
            "FokI" => "AGIJKMNQRVWX",
            "FriOI" => "IV",
            "FseI" => "AN",
            "FspI" => "JNO",
            "FspAI" => "F",
            "FspBI" => "F",
            "Fsp4HI" => "I",
            "GlaI" => "I",
            "GluI" => "I",
            "GsuI" => "F",
            "HaeII" => "GJKMNORSW",
            "HaeIII" => "ABGHIJKMNOQRSUWXY",
            "HapII" => "AK",
            "HgaI" => "IN",
            "HhaI" => "ABFGJKNORUWY",
            "Hin1I" => "FKO",
            "Hin1II" => "F",
            "Hin4I" => "F",
            "Hin6I" => "F",
            "HinP1I" => "N",
            "HincII" => "ABFGHJKNOQRUWXY",
            "HindII" => "IMSV",
            "HindIII" => "ABCFGHIJKMNOQRSUVWXY",
            "HinfI" => "ABCFGHIJKMNOQRUVWXY",
            "HpaI" => "ABCGHIJKMNOQRSUVWX",
            "HpaII" => "BFGIMNOQRSUVWX",
            "HphI" => "FN",
            "Hpy8I" => "F",
            "Hpy99I" => "N",
            "Hpy188I" => "N",
            "Hpy188III" => "N",
            "HpyAV" => "N",
            "HpyCH4III" => "N",
            "HpyCH4IV" => "N",
            "HpyCH4V" => "N",
            "HpyF3I" => "F",
            "HpyF10VI" => "F",
            "Hsp92I" => "R",
            "Hsp92II" => "R",
            "HspAI" => "IV",
            "ItaI" => "M",
            "KasI" => "N",
            "KpnI" => "ABCFGHIJKMNOQRSUVWXY",
            "Kpn2I" => "F",
            "KspI" => "MS",
            "Ksp22I" => "IV",
            "Ksp632I" => "M",
            "KspAI" => "F",
            "Kzo9I" => "I",
            "LguI" => "F",
            "LweI" => "F",
            "MabI" => "I",
            "MaeI" => "M",
            "MaeII" => "M",
            "MaeIII" => "M",
            "MalI" => "I",
            "MamI" => "M",
            "MbiI" => "F",
            "MboI" => "ABCFGKNQRUWXY",
            "MboII" => "AFGIJKNOQRVWX",
            "MfeI" => "N",
            "MflI" => "K",
            "MhlI" => "IV",
            "MlsI" => "F",
            "MluI" => "ABFGHIJKMNOQRSUVWX",
            "MluNI" => "MS",
            "MlyI" => "N",
            "Mly113I" => "I",
            "MmeI" => "NX",
            "MnlI" => "FGINQVWX",
            "Mph1103I" => "F",
            "MreI" => "F",
            "MroI" => "MO",
            "MroNI" => "IV",
            "MroXI" => "IV",
            "MscI" => "BNO",
            "MseI" => "BN",
            "MslI" => "N",
            "MspI" => "AFGHIJKMNOQRSUVWXY",
            "Msp20I" => "IV",
            "MspA1I" => "INRV",
            "MspCI" => "C",
            "MspR9I" => "I",
            "MssI" => "F",
            "MunI" => "FKM",
            "MvaI" => "AFGKMOSW",
            "Mva1269I" => "F",
            "MvnI" => "M",
            "MvrI" => "U",
            "MwoI" => "N",
            "NaeI" => "ACKMNORU",
            "NarI" => "GJMNOQRUWX",
            "NciI" => "GJNORSW",
            "NcoI" => "ABCFGHJKMNOQRSUWXY",
            "NdeI" => "ABFGJKMNQRSWXY",
            "NdeII" => "GJMRSW",
            "NgoMIV" => "NR",
            "NheI" => "ABFGJKMNORSUW",
            "NlaIII" => "GNW",
            "NlaIV" => "GNW",
            "NmeAIII" => "N",
            "NmuCI" => "F",
            "NotI" => "ABCFGHJKMNOQRSUWXY",
            "NruI" => "ABCGIJKMNOQRSUWX",
            "NsbI" => "FK",
            "NsiI" => "BGHJMNRSUW",
            "NspI" => "MN",
            "NspV" => "JO",
            "OliI" => "F",
            "PacI" => "GNOW",
            "PaeI" => "F",
            "PaeR7I" => "N",
            "PagI" => "F",
            "PalAI" => "I",
            "PasI" => "F",
            "PauI" => "F",
            "PceI" => "IV",
            "PciI" => "IN",
            "PciSI" => "I",
            "PctI" => "IV",
            "PdiI" => "F",
            "PdmI" => "F",
            "PfeI" => "F",
            "Pfl23II" => "F",
            "PflFI" => "N",
            "PflMI" => "N",
            "PfoI" => "F",
            "PhoI" => "N",
            "PinAI" => "BM",
            "PleI" => "N",
            "Ple19I" => "I",
            "PmaCI" => "AK",
            "PmeI" => "GNW",
            "PmlI" => "N",
            "PpiI" => "F",
            "PpsI" => "I",
            "Ppu21I" => "F",
            "PpuMI" => "NO",
            "PscI" => "F",
            "PshAI" => "AKN",
            "PshBI" => "K",
            "PsiI" => "IN",
            "Psp5II" => "F",
            "Psp6I" => "I",
            "Psp1406I" => "FK",
            "Psp124BI" => "IV",
            "PspCI" => "IV",
            "PspEI" => "IV",
            "PspGI" => "N",
            "PspLI" => "I",
            "PspN4I" => "I",
            "PspOMI" => "INV",
            "PspPPI" => "I",
            "PspXI" => "IN",
            "PsrI" => "I",
            "PstI" => "ABCFGHIJKMNOQRSUVWXY",
            "PsuI" => "F",
            "PsyI" => "F",
            "PvuI" => "ABFGKMNOQRSUWXY",
            "PvuII" => "ABCFGHIJKMNORSUVWXY",
            "RcaI" => "M",
            "RgaI" => "I",
            "RigI" => "I",
            "RsaI" => "BCFGHIJMNOQRSVWXY",
            "RsaNI" => "I",
            "RseI" => "F",
            "RsrII" => "MNQX",
            "Rsr2I" => "IV",
            "SacI" => "AFGHJKMNOQRSUWX",
            "SacII" => "AGHJKNOQRWX",
            "SalI" => "ABCFGHIJKMNOQRSUVWXY",
            "SanDI" => "E",
            "SapI" => "N",
            "SatI" => "F",
            "Sau96I" => "GJMNOUW",
            "Sau3AI" => "AGHJKMNOQRSUWX",
            "SbfI" => "INV",
            "ScaI" => "ABCFGJKMNOQRSWX",
            "SchI" => "F",
            "ScrFI" => "JMNOS",
            "SdaI" => "F",
            "SduI" => "F",
            "SetI" => "I",
            "SexAI" => "MN",
            "SfaNI" => "INV",
            "SfcI" => "N",
            "SfiI" => "ACFGIJKMNOQRSUVWX",
            "SfoI" => "N",
            "Sfr274I" => "IV",
            "Sfr303I" => "IV",
            "SfuI" => "M",
            "SgfI" => "R",
            "SgrAI" => "MN",
            "SgrBI" => "C",
            "SgrDI" => "F",
            "SgsI" => "F",
            "SinI" => "GQRWX",
            "SlaI" => "C",
            "SmaI" => "ABCFGHIJKMNOQRSUVWXY",
            "SmiI" => "FIKV",
            "SmiMI" => "IV",
            "SmlI" => "N",
            "SmoI" => "F",
            "SmuI" => "F",
            "SnaBI" => "ACKMNR",
            "SpeI" => "ABGHJKMNOQRSUWX",
            "SphI" => "ABCGHIJKMNOQRSVWX",
            "SrfI" => "EO",
            "Sse9I" => "IV",
            "Sse8387I" => "AK",
            "SseBI" => "C",
            "SsiI" => "F",
            "SspI" => "ABCFGIJKMNOQRSUVWX",
            "SstI" => "BC",
            "SstII" => "B",
            "StrI" => "U",
            "StuI" => "ABJKMNQRSUX",
            "StyI" => "CJMNRS",
            "StyD4I" => "N",
            "SwaI" => "GJMNSW",
            "TaaI" => "F",
            "TaiI" => "F",
            "TaqI" => "ABCFGIJKMNOQRSUVWXY",
            "TaqII" => "QX",
            "TasI" => "F",
            "TatI" => "F",
            "TauI" => "F",
            "TfiI" => "N",
            "TliI" => "N",
            "Tru1I" => "F",
            "Tru9I" => "GIMRVW",
            "TseI" => "N",
            "TsoI" => "F",
            "Tsp45I" => "N",
            "Tsp509I" => "N",
            "TspDTI" => "X",
            "TspEI" => "O",
            "TspGWI" => "X",
            "TspMI" => "N",
            "TspRI" => "N",
            "TstI" => "F",
            "Tth111I" => "GIKNQRVWX",
            "Van91I" => "AFKM",
            "Vha464I" => "IV",
            "VneI" => "IV",
            "VpaK11BI" => "K",
            "VspI" => "FIRV",
            "XagI" => "F",
            "XapI" => "F",
            "XbaI" => "ABCFGHIJKMNOQRSUVWXY",
            "XceI" => "F",
            "XcmI" => "N",
            "XhoI" => "ABFGHJKMNOQRSUWXY",
            "XhoII" => "GMRW",
            "XmaI" => "INRUV",
            "XmaCI" => "M",
            "XmaJI" => "F",
            "XmiI" => "F",
            "XmnI" => "GNRUW",
            "XspI" => "K",
            "ZraI" => "INV",
            "ZrmI" => "I",
            "Zsp2I" => "IV"
        );
        return $vendors;
    }

    private function extractSequence($text) {
        $sequence["seq"] = preg_replace("/\W|\d/", "", strtoupper($text));
        return $sequence;
    }
    
}

?>