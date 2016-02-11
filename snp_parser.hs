module Snp_parser where

import Snp_windows

import Control.Monad.Trans.Except
import qualified Data.Text.Lazy as T
import qualified Data.Text.Lazy.IO as TIO
import Data.Text.Lazy.Read

data ParseError = UnsortedSNPs 
                | PositionUnparsed String
                | ValueUnparsed String
                | TooFewFields Int Int
                | ParseError String
                deriving (Show, Eq)



type ParsedSNP = Except ParseError (SNP Double)

parseSNPlist :: [T.Text] -> [ParsedSNP]
parseSNPlist [] = []
parseSNPlist (l:[]) = [parse l]
parseSNPlist (l:lst) = let
  parsed1 = runExcept $ parse l
  parsed2 = runExcept . parse $ head lst
  parsed  = except $ isOrdered parsed1 parsed2 -- if both p1 & p2 are valid: check order
  in parsed : parseSNPlist lst
    where
      isOrdered (Left x) _          = Left x
      isOrdered (Right x) (Right y) 
           | pos x > pos y  &&
             chrom x == chrom y     = Left UnsortedSNPs
           | otherwise              = Right x
      isOrdered (Right x) _         = Right x

{- Main parser function,
   for expanding to multiple values, edit this function -}
parse :: T.Text -> ParsedSNP
parse txt = let
  fields = T.splitOn (T.pack "\t") txt
  in if not $ fields `longerThan` 2
     then except . Left $ TooFewFields (length fields) 3
     else toSNP fields

  where
    toSNP (x1:x2:x3:xs) = do
      pos <- catchE (except $ decimal x2) 
             (\e -> throwE . PositionUnparsed $ "error parsing position-field: " ++ e)
      val <- catchE (except $ double x3) 
             (\e -> throwE . ValueUnparsed $ "error parsing value-field: " ++ e)
      return $ SNP (T.unpack x1) (fst pos) (fst val)

longerThan :: [a] -> Int -> Bool
longerThan lst len = let
  mlen = dropWhile ( <= len) runningLength
  runningLength = (map snd $ zip lst [1..])
  in case mlen of
      [] -> False
      xs -> head xs > len
