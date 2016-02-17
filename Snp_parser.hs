module Snp_parser (parseSNPs, ParseError (..), reportSNPerror) where

import Snp_windows

import Data.Maybe                    (catMaybes)
import System.IO                     (hPutStrLn, stderr)
import System.Exit                   (exitFailure)
import Control.Monad.Trans.Except    (Except, except, runExcept, catchE, throwE)
import Data.Text.Lazy                (Text, pack, unpack, splitOn)
import Data.Text.Lazy.Read           (decimal, double)

data ParseError = UnsortedSNPs 
                | PositionUnparsed String
                | ValueUnparsed String
                | TooFewFields Int Int
                | ParseError String
                deriving (Show, Eq)

type ParsedSNP a = Except ParseError (SNP a)

parseSNPs :: [Text] -> [ParsedSNP Double]
parseSNPs = assertSNPorder . (map parseFields)

reportSNPs :: [ParsedSNP a] -> IO [SNP a]
reportSNPs parsed_snps = do
  m_snps <- sequence . (map reportSNPerror) $ enumerate parsed_snps
  return $ catMaybes m_snps
  where
    enumerate = zip [1..]

-- take line-number and parsed SNP and performs IO to report parse-error to user.
reportSNPerror :: (Int, ParsedSNP a) -> IO (Maybe (SNP a))
reportSNPerror (i,x) = 
  case runExcept x of
    Left UnsortedSNPs -> do 
      printErr $ "The SNPs are not in order (line " ++ show i ++
                  " and " ++ show (i + 1) ++ ")."
    Left (PositionUnparsed e) -> do
      printErr $ "The position-field could not be parsed at line " ++ show i ++ ": " ++ e
    Left (ValueUnparsed e) -> do
      printErr $ "Could not parse the value-field(s) at line " ++ show i ++": " ++ e
    Left (TooFewFields a b) -> do
      printErr $ "Missing fields at line " ++ show i ++ ". " ++ show a ++ " fields present when " ++
                     show b ++ " fields is the minimum required"
    Left (ParseError e) -> do
      printErr $ "Couldn't parse line " ++ show i ++ ": " ++ e
    Right x -> return $ Just x
        
  where
    printErr e = do 
      hPutStrLn stderr $ "Error: " ++ e
      exitFailure
      return Nothing

assertSNPorder :: [ParsedSNP a] -> [ParsedSNP a]
assertSNPorder = (map except) . (foldr f []) . (map runExcept)
  where
    f (Right x) l@((Right y):_) 
      | pos x > pos y && chrom x == chrom y = (Left UnsortedSNPs) : l
      | otherwise                           = (Right x) : l
    f x l = x : l

{-
   Main parser function.
   For expanding to multiple values, edit this function
   along with the ParsedSNP-alias 
 -}
parseFields :: Text -> ParsedSNP Double
parseFields txt = let
  fields = splitOn (pack "\t") txt
  in if not $ fields `longerThan` 2
     then except . Left $ TooFewFields (length fields) 3
     else toSNP fields
  where
    toSNP (x1:x2:x3:_) = do
      pos <- catchE (except $ decimal x2) (\e -> throwE $ PositionUnparsed e)
      val <- catchE (except $ double x3) (\e -> throwE $ ValueUnparsed e)
      return $ SNP x1 (fst pos) (fst val)

longerThan :: [a] -> Int -> Bool
longerThan lst len = let
  mlen = dropWhile ( <= len) runningLength
  runningLength = (map snd $ zip lst [1..])
  in case mlen of
      [] -> False
      xs -> head xs > len

