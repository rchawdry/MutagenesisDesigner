<?php

    header('Content-Type: application/json');
    
    $dna = isset($_REQUEST['dna']) ? $_REQUEST['dna']:'';
    $mutationPos = isset($_REQUEST['mutationPos']) ? $_REQUEST['mutationPos']:'';
    $mutationLen = isset($_REQUEST['mutationLen']) ? $_REQUEST['mutationLen']:'';
    $numPrimers = isset($_REQUEST['numPrimers']) ? $_REQUEST['numPrimers']:'1';
    $action = isset($_REQUEST['action']) ? $_REQUEST['action']:'';
    $returnArray = isset($_REQUEST['returnArray']) ? $_REQUEST['returnArray']:'';    

    $baseUrl="http://localhost:81/webservice/";
    //$baseUrl="http://www.vixis.com/webservice/";
    $enzymeUrl = $baseUrl."enzymedesign.php";
    $primerUrl = $baseUrl."primerdesign.php";
    

    $result=array();
    
    if($action=='' || $action='e'){
        $enzymeList = getResult($enzymeUrl, $dna, $mutationPos, $mutationLen, $numPrimers, $returnArray);
        $result['enzymes']=json_decode($enzymeList,true);
    }
    if($action=='' || $action='p'){
        $primerList = getResult($primerUrl, $dna, $mutationPos, $mutationLen, $numPrimers,  $returnArray);
        $result['primers']=json_decode($primerList,true);
    }

    echo json_encode($result);

    function getResult($url, $dna, $mutationPos, $mutationLen, $numPrimers, $returnArray){
        $ch = curl_init($url);
        $options = array(
            CURLOPT_RETURNTRANSFER => true,
            CURLOPT_POST => 1,
            CURLOPT_POSTFIELDS => array(
                'dna'=>$dna,
                'mutationPos'=>$mutationPos,
                'mutationLen'=>$mutationLen,
                'numPrimers'=>$numPrimers,
                'returnArray'=>$returnArray,
                'jsonoutput'=>1,
                )
        );
        curl_setopt_array($ch, $options);
        return curl_exec($ch); 
    }
    
?>
