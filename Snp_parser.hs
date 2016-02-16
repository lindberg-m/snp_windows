module Snp_parser (parseSNPlist, ParseError (..)) where

import Snp_windows

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

parseSNPlist :: [Text] -> [ParsedSNP Double]
parseSNPlist = assertSNPorder . (map parse)

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
parse :: Text -> ParsedSNP Double
parse txt = let
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

{-
  Unused function, for general parsing instead of "parseSNPlist" 
  if multiple values are present
-}
parseDoubles :: [Text] -> [ParsedSNP [Double]]
parseDoubles = assertSNPorder . map (parse' $ fmap fst . double)

parse' :: (Text -> Either String a) -> Text -> ParsedSNP [a]
parse' f txt = let
  fields = splitOn (pack "\t") txt
  in if not $ fields `longerThan` 2
     then except . Left $ TooFewFields (length fields) 3
     else toSNP fields
  where
    toSNP (x1:x2:xs) = do
      pos <- catchE (except $ decimal x2) (\e -> throwE $ PositionUnparsed e)
      vals <- parseValues f xs 
      return $ SNP x1 (fst pos) vals

    parseValues :: (Text -> Either String a) -> [Text] -> Except ParseError [a]
    parseValues f = sequence . map (\x -> catchE (except $ f x) (\e -> throwE $ ValueUnparsed e))


