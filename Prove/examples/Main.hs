{-# LANGUAGE NamedFieldPuns #-}
module Main where

import Protolude
-- import Control.Monad.Random (getRandomR)
import Bulletproofs.ArithmeticCircuit
import Data.Pairing.BN254 (Fr)
-- import Data.Field.Galois (rnd)
import Data.List.Split (divvy)

-- import Sonic.SRS as SRS
import Sonic.Protocol
import Sonic.Circuits
-- import Sonic.CommitmentScheme (pcVShow)
-- import Sonic.Signature (hscVerifyShow)
-- import Sonic.Constraints (sPoly)

import qualified Data.Text    as Text
import qualified Data.Text.IO as Text
import Data.Text.Read (decimal, signed)

import Data.Time
-- import Control.DeepSeq (force)

-- import qualified Data.ByteString
-- import qualified Crypto.Hash.SHA256 as SHA256

-- sonicProtocol :: SRS 
--   -> SRS 
--   -> Proof 
--   -> RndOracle 
--   -> ArithCircuit Fr 
--   -> Assignment Fr 
--   -> Fr 
--   -> IO Bool
-- sonicProtocol srsRaw srsLocal proof rndOracle@RndOracle{..} circuit assignment x = do
--   -- Verifier
--   pure $ verify srsRaw srsLocal circuit proof rndOracleY rndOracleZ rndOracleYZs


outputProof:: ArithCircuit Fr -> Assignment Fr -> Fr -> Fr -> Fr -> Fr -> IO ()
outputProof circuit assignment pXRaw pXLocal alphaRaw alphaLocal = do
  -- Setup for an SRS
  -- startSrs <- getCurrentTime
  -- srsRaw <- SRS.new <$> pure n <*> pure pXRaw <*> pure alphaRaw
  -- srsLocal <- SRS.new <$> pure n <*> pure pXLocal <*> pure alphaLocal

  -- print $ alphaRaw
  -- print $ pXRaw
  -- print $ alphaLocal
  -- print $ (srsD srsRaw)
  -- print $ (srsD srsLocal)
  -- print $ "writing SRSs"
  -- writeFile "output/srsRaw.txt" $ show $ srsRaw
  -- writeFile "output/srsLocal.txt" $ show $ srsLocal
  -- stopSrs <- getCurrentTime
  -- print $ diffUTCTime stopSrs startSrs
  -- print $ pXRaw
  -- print $ pXLocal
  -- Prover
  -- (proofOutsourced, rndOracle1) <- proveOutsourced srsLocal srsLocal 4 assignment circuit
  -- putText $ "proofOutsourced: " <> show proofOutsourced
  -- print $ "generating proof:"
  start <- getCurrentTime
  (proof) <- prove 4 n assignment circuit
  -- putText $ "proof: " <> show proof
  -- putText $ "polys: " <> show verifierData
  writeFile "output/polys_dt.txt" $ show $ proof
  stop <- getCurrentTime
  print $ diffUTCTime stop start
  -- print $ "verifying proof:"
  -- startVer <- getCurrentTime
  -- putText $ "success:" <> show (verify srsLocal srsLocal circuit proof proofOutsourced rndOracleY rndOracleZ rndOracleYZ)
  -- stopVer <- getCurrentTime
  -- print $ diffUTCTime stopVer startVer

  -- putText $ show (hscVerifyShow srsLocal sXY rndOracleYZs prHscProof)

  -- putText $ show (pcVShow srsRaw (fromIntegral nexample) prRRaw rndOracleZ (prARaw, prWaRaw))
  -- print $ "writing verify Hxi String"
  -- writeFile "output/verifyHxiString.txt" $ show (getverifyHxiString srsRaw srsLocal circuit proof rndOracleY rndOracleZ rndOracleYZs)
  -- print $ "writing proof"
  -- writeFile "output/proof.txt" $ show $ proof
  -- print $ "writing rndOracle"
  -- writeFile "output/rndOracle.txt" $ show $ rndOracle
  -- print $ "writing verifier data"
  -- writeFile "output/verifierData.txt" $ show $ verifierData
  where
    -- n: Number of multiplication constraints
    n = (length $ aL assignment) * 6
    -- nexample = 50
    -- sXY = sPoly (weights circuit)
    


runExample :: IO ()
runExample = do

-- alphaRaw = 4537460542209314651160888417413866091249215769242027952878258319870902529429
--       alphaLocal = 10318208573435976958553514324035131570702768863559211839926632048062238413316
--       pXRaw = 4518563069472097478295524977775021906947577384653869551543466909390271555451
--       pXLocal = 19708723214916757413126173169122466312825114221484651297201668130676937834219



  let alphaRaw = 10
      alphaLocal = 11
      pXRaw = 12
      pXLocal = 13

  wLS <- fmap Text.words (Text.readFile "input/wL.txt")
  wRS <- fmap Text.words (Text.readFile "input/wR.txt")
  wOS <- fmap Text.words (Text.readFile "input/wO.txt")
  csS <- fmap Text.words (Text.readFile "input/cs.txt")
  aLS <- fmap Text.words (Text.readFile "input/aL.txt")
  aRS <- fmap Text.words (Text.readFile "input/aR.txt")
  aOS <- fmap Text.words (Text.readFile "input/aO.txt")
  let wLL = foldr (\x acc -> (fst x):acc) [] (rights (map (signed decimal) wLS))
      wRL = foldr (\x acc -> (fst x):acc) [] (rights (map (signed decimal) wRS))
      wOL = foldr (\x acc -> (fst x):acc) [] (rights (map (signed decimal) wOS))
      cs = foldr (\x acc -> (fst x):acc) [] (rights (map (signed decimal) csS))
      aL = foldr (\x acc -> (fst x):acc) [] (rights (map (signed decimal) aLS))
      aR = foldr (\x acc -> (fst x):acc) [] (rights (map (signed decimal) aRS))
      aO = foldr (\x acc -> (fst x):acc) [] (rights (map (signed decimal) aOS))
      inputSize = 232
      -- wL = divvy 50 50 wLL
      -- wR = divvy 50 50 wRL
      -- wO = divvy 50 50 wOL
      wL = divvy inputSize inputSize wLL
      wR = divvy inputSize inputSize wRL
      wO = divvy inputSize inputSize wOL

      (arithCircuit, assignment) = arithCircuitExample wL wR wO cs aL aR aO
  -- success <- sonicProtocol arithCircuit assignment pX
  -- putText $ "Success: " <> show success
  outputProof arithCircuit assignment pXRaw pXLocal alphaRaw alphaLocal



main :: IO ()
main = runExample
-- main = do
--   wLS <- fmap Text.words (Text.readFile "input/sample_wL.txt")
--   wRS <- fmap Text.words (Text.readFile "input/sample_wR.txt")
--   wOS <- fmap Text.words (Text.readFile "input/sample_wO.txt")
--   csS <- fmap Text.words (Text.readFile "input/sample_cs.txt")
--   aLS <- fmap Text.words (Text.readFile "input/sample_aL.txt")
--   aRS <- fmap Text.words (Text.readFile "input/sample_aR.txt")
--   aOS <- fmap Text.words (Text.readFile "input/sample_aO.txt")
--   let wLL = foldr (\x acc -> (fst x):acc) [] (rights (map (signed decimal) wLS))
--       wRL = foldr (\x acc -> (fst x):acc) [] (rights (map (signed decimal) wRS))
--       wOL = foldr (\x acc -> (fst x):acc) [] (rights (map (signed decimal) wOS))
--       cs = foldr (\x acc -> (fst x):acc) [] (rights (map (signed decimal) csS))
--       aL = foldr (\x acc -> (fst x):acc) [] (rights (map (signed decimal) aLS))
--       aR = foldr (\x acc -> (fst x):acc) [] (rights (map (signed decimal) aRS))
--       aO = foldr (\x acc -> (fst x):acc) [] (rights (map (signed decimal) aOS))
--       wL = divvy 2 2 wLL
--       wR = divvy 2 2 wRL
--       wO = divvy 2 2 wOL
--   print wL
--   print wR
--   print wO
--   print cs
--   print aL
--   print aR
--   print aO
  

