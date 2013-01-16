<?php
    header('Content-Type: application/json');

    $dna = isset($_REQUEST['dna']) ? $_REQUEST['dna']:'';
    $mutationPos = isset($_REQUEST['mutationPos']) ? $_REQUEST['mutationPos']:'';
    $mutationLen = isset($_REQUEST['mutationLen']) ? $_REQUEST['mutationLen']:'';
    $numPrimers = isset($_REQUEST['numPrimers']) ? $_REQUEST['numPrimers']:5;
    $jsonOutputFlg = isset($_REQUEST['jsonoutput']) ? 1:0;	
    
    if (strlen($dna)>0)
    {
        $primer_output = designPrimer($dna, $mutationPos,$mutationLen, $numPrimers);
        $primer_array = explodeBoulderIO($primer_output);
        if (count($primer_array)>0)
        { 
            if ($jsonOutputFlg)
                echo json_encode(getStructuredArray($primer_array));
            else
                echo json_encode($primer_array);
        }
        else
        {
            echo $primer_output;   
        }
    }
    else
    {
        echo "Please submit a DNA Sequence using the 'dna' form parameter";

    }

    
function designPrimer($dnaSeq, $mutationPos, $mutationLen, $numPrimers) {

    //Need to create a unique session identifier to create temp files for primer software
    $session_id = rand_string(15);

    //Constants for location of Primer3 software

    //Windows
    $swPath = 'C:/Sypherlink/Rehan/Vixis/FirstBioSciences/Software/primer3/primer3-2.3.4/test2/';
    $sessionPath = $swPath . 'sessions/';
    $swCmd = $swPath . 'primer3_core.exe  -p3_settings_file=' . $swPath . 'settings.txt';


    //Linux
    /*
    $swPath = './primer3/2.3.4/';
    $sessionPath = './sessions/';
    $swCmd = $swPath . 'primer3_core -p3_settings_file=' . $swPath . 'settings.txt';
     */
     
    $params = array();

    //Fill out with passed-in parameters
    $params['SEQUENCE_ID']= $session_id;
    $params['SEQUENCE_TEMPLATE']=$dnaSeq;
    $params['PRIMER_PICK_RIGHT_PRIMER']=1;
    $params['PRIMER_PICK_INTERNAL_OLIGO']=0;
    $params['PRIMER_PICK_LEFT_PRIMER']=1;
    $params['PRIMER_NUM_RETURN']=$numPrimers;
    $params['PRIMER_TASK']='pick_detection_primers';
    $params['PRIMER_PRODUCT_SIZE_RANGE']='150-300';

    if($mutationPos!=''){
        $params["SEQUENCE_TARGET"]=$mutationPos.','.$mutationLen;
    }
    
    //Save the file to be used as input to the primer software
    $newFileText=implodeBoulderIO($params);
    $swInputFile = $sessionPath . "input_" . $params["SEQUENCE_ID"] . ".txt";
    file_put_contents($swInputFile, $newFileText);
    
    //Execute shell script
    $swExe = $swCmd . ' < ' . $swInputFile . ' 2>&1';
    $output = shell_exec($swExe);
 
    unlink($swInputFile);
		
    return $output;
}

function explodeBoulderIO($boulderIOString){
    $pairs=array();
    preg_match_all("/([^=]+)=(.*)\n/",$boulderIOString,$pairs);
    $content_array = array_combine($pairs[1], $pairs[2]);
    return $content_array;
}

function implodeBoulderIO($input_array){
    $content_string='';
    foreach ($input_array as $key => $value) {
        $content_string .= $key."=".$value."\n";
    }
    return $content_string.="=\n";
}
function rand_string( $length ) {
	$chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
        $str="";
	$size = strlen( $chars );
	for( $i = 0; $i < $length; $i++ ) {
		$str .= $chars[ rand( 0, $size - 1 ) ];
	}

	return $str;
}
function getStructuredArray( $input ) {
	$o=array(
            'input'=>array(
                'SEQUENCE_TEMPLATE'=>$input['SEQUENCE_TEMPLATE'],
                'PRIMER_PICK_RIGHT_PRIMER'=>$input['PRIMER_PICK_RIGHT_PRIMER'],
                'PRIMER_PICK_INTERNAL_OLIGO'=>$input['PRIMER_PICK_INTERNAL_OLIGO'],
                'PRIMER_PICK_LEFT_PRIMER'=>$input['PRIMER_PICK_LEFT_PRIMER'],
                'PRIMER_NUM_RETURN'=>$input['PRIMER_NUM_RETURN'],
                'PRIMER_TASK'=>$input['PRIMER_TASK'],
                ),
            'output'=>array(
                'PRIMER_LEFT_NUM_RETURNED'=>$input['PRIMER_LEFT_NUM_RETURNED'],
                'PRIMER_RIGHT_NUM_RETURNED'=>$input['PRIMER_RIGHT_NUM_RETURNED'],
                'PRIMER_INTERNAL_NUM_RETURNED'=>$input['PRIMER_INTERNAL_NUM_RETURNED'],
                'PRIMER_PAIR_NUM_RETURNED'=>$input['PRIMER_PAIR_NUM_RETURNED'],
            ),
        );
        
        $primers=array();
        $numPrimers = $input['PRIMER_LEFT_NUM_RETURNED'];
        for($i=0;$i<$numPrimers;$i++){
            $primerLeftPosLen = explode(',',$input['PRIMER_LEFT_' . $i]);
            $primerRightPosLen = explode(',',$input['PRIMER_RIGHT_' . $i]);
            
            $primer = array(
                'PRIMER_PAIR_PENALTY'=>$input['PRIMER_PAIR_' . $i. '_PENALTY'],
                'PRIMER_LEFT_PENALTY'=>$input['PRIMER_LEFT_' . $i. '_PENALTY'],
                'PRIMER_RIGHT_PENALTY'=>$input['PRIMER_RIGHT_' . $i. '_PENALTY'],
                'PRIMER_LEFT_SEQUENCE'=>$input['PRIMER_LEFT_' . $i. '_SEQUENCE'],
                'PRIMER_RIGHT_SEQUENCE'=>$input['PRIMER_RIGHT_' . $i. '_SEQUENCE'],
                'PRIMER_LEFT'=>$input['PRIMER_LEFT_' . $i],
                'PRIMER_RIGHT'=>$input['PRIMER_RIGHT_' . $i],
                'PRIMER_LEFT_TM'=>$input['PRIMER_LEFT_' . $i. '_TM'],
                'PRIMER_RIGHT_TM'=>$input['PRIMER_RIGHT_' . $i. '_TM'],
                'PRIMER_LEFT_GC_PERCENT'=>$input['PRIMER_LEFT_' . $i. '_GC_PERCENT'],
                'PRIMER_RIGHT_GC_PERCENT'=>$input['PRIMER_RIGHT_' . $i. '_GC_PERCENT'],
                'PRIMER_LEFT_SELF_ANY_TH'=>$input['PRIMER_LEFT_' . $i. '_SELF_ANY_TH'],
                'PRIMER_RIGHT_SELF_ANY_TH'=>$input['PRIMER_RIGHT_' . $i. '_SELF_ANY_TH'],
                'PRIMER_LEFT_SELF_END_TH'=>$input['PRIMER_LEFT_' . $i. '_SELF_END_TH'],
                'PRIMER_RIGHT_SELF_END_TH'=>$input['PRIMER_RIGHT_' . $i. '_SELF_END_TH'],
                'PRIMER_LEFT_HAIRPIN_TH'=>$input['PRIMER_LEFT_' . $i. '_HAIRPIN_TH'],
                'PRIMER_RIGHT_HAIRPIN_TH'=>$input['PRIMER_RIGHT_' . $i. '_HAIRPIN_TH'],
                'PRIMER_LEFT_END_STABILITY'=>$input['PRIMER_LEFT_' . $i. '_END_STABILITY'],
                'PRIMER_RIGHT_END_STABILITY'=>$input['PRIMER_RIGHT_' . $i. '_END_STABILITY'],
                'PRIMER_PAIR_COMPL_ANY_TH'=>$input['PRIMER_PAIR_' . $i. '_COMPL_ANY_TH'],
                'PRIMER_PAIR_COMPL_END_TH'=>$input['PRIMER_PAIR_' . $i. '_COMPL_END_TH'],
                'PRIMER_PAIR_PRODUCT_SIZE'=>$input['PRIMER_PAIR_' . $i. '_PRODUCT_SIZE'],
                'PRIMER_LEFT_POSITION'=>$primerLeftPosLen[0],
                'PRIMER_LEFT_LENGTH'=>$primerLeftPosLen[1],
                'PRIMER_RIGHT_POSITION'=>$primerRightPosLen[0],
                'PRIMER_RIGHT_LENGTH'=>$primerRightPosLen[1],
            );
            array_push($primers, $primer);
        }
        $o['primers']=$primers;
	return $o;
}
?>
