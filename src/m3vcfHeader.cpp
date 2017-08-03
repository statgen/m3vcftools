#include "m3vcfHeader.h"

m3vcfHeader::m3vcfHeader() : myHeaderLines()
{
    reset();
}


m3vcfHeader::~m3vcfHeader()
{
}


//bool m3vcfHeader::read(IFILE filePtr)
//{
//    // Reading, so clean out this header.
//    reset();
//
//    if(filePtr == NULL)
//    {
//        // No file was passed in.
//        myStatus.setStatus(StatGenStatus::FAIL_ORDER,
//                           "Need to pass in an open file ptr to m3vcfHeader::read.");
//        return(false);
//    }
//
//    // Read until the header line has been read (after the meta lines).
//    while(!myHasHeaderLine)
//    {
//        // Increase the size of headerlines by 1 to fit the new line.
//        myHeaderLines.resize(myHeaderLines.size() + 1);
//
//        // Read the next line from the file into the header structure.
//        String& newStr = myHeaderLines.back();
//        if(newStr.ReadLine(filePtr) < 0)
//        {
//            // Error, unable to read an entire header from the file.
//            myStatus.setStatus(StatGenStatus::INVALID,
//                               "Error reading M3VCF Meta/Header, EOF found before the header line.");
//            return(false);
//        }
//        if(newStr.Length() <= 2)
//        {
//            // A header/meta line must have at least 2 characters
//            // ## or # and 8 fields, so if less than 2 characters,
//            // error.
//            myStatus.setStatus(StatGenStatus::INVALID,
//                               "Error reading VCF Meta/Header, line without at least 2 characters found before the header line.");
//            return(false);
//        }
//
//        // Check if it is a header (first char is # and 2nd one is not).
//        if((newStr[0] == '#') && (newStr[1] != '#'))
//        {
//            myHasHeaderLine = true;
//
//            // Parse the header line to get the sample information.
//            myParsedHeaderLine.ReplaceColumns(newStr, '\t');
//            return SparseSampleNames();
//        }
//        else if((newStr[0] != '#') || (newStr[1] != '#'))
//        {
//            // A meta line must start with "##", we expect meta lines until
//            // the header line is found.
//            myStatus.setStatus(StatGenStatus::INVALID,
//                               "Error reading VCF Meta/Header, line not starting with '##' found before the header line.");
//            return(false);
//        }
//    }
//    return(true);
//}

bool m3vcfHeader::read(IFILE filePtr)
{
    // Reading, so clean out this header.
    reset();

    if(filePtr == NULL)
    {
        // No file was passed in.
        myStatus.setStatus(StatGenStatus::FAIL_ORDER,
                           "Need to pass in an open file ptr to m3vcfHeader::read.");
        return(false);
    }

    // Read until the header line has been read (after the meta lines).
    while(!myHasHeaderLine)
    {
        // Increase the size of headerlines by 1 to fit the new line.
        myHeaderLinesNew.resize(myHeaderLinesNew.size() + 1);

        // Read the next line from the file into the header structure.
        string& newStr = myHeaderLinesNew.back();
        if(filePtr->readLine(newStr) < 0)
        {
            // Error, unable to read an entire header from the file.
            myStatus.setStatus(StatGenStatus::INVALID,
                               "Error reading M3VCF Meta/Header, EOF found before the header line.");
            return(false);
        }
        if(newStr.length() <= 2)
        {
            // A header/meta line must have at least 2 characters
            // ## or # and 8 fields, so if less than 2 characters,
            // error.
            myStatus.setStatus(StatGenStatus::INVALID,
                               "Error reading VCF Meta/Header, line without at least 2 characters found before the header line.");
            return(false);
        }

        // Check if it is a header (first char is # and 2nd one is not).
        if((newStr[0] == '#') && (newStr[1] != '#'))
        {
            myHasHeaderLine = true;

            // Parse the header line to get the sample information.
            String tempString = newStr.c_str();
            myParsedHeaderLine.ReplaceColumns(tempString, '\t');
            return SparseSampleNames();
        }
        else if((newStr[0] != '#') || (newStr[1] != '#'))
        {
            // A meta line must start with "##", we expect meta lines until
            // the header line is found.
            myStatus.setStatus(StatGenStatus::INVALID,
                               "Error reading VCF Meta/Header, line not starting with '##' found before the header line.");
            return(false);
        }
    }
    return(true);
}


void m3vcfHeader::reset()
{
    myHasHeaderLine = false;
    myHeaderLines.clear();
    myHeaderLinesNew.clear();
    SampleNames.clear();
    SamplePloidy.clear();
    numSamples=0;
}


bool m3vcfHeader::SparseSampleNames()
{
    int position = NUM_M3VCF_NON_SAMPLE_HEADER_COLS;

    while(position<myParsedHeaderLine.Length())
    {
        SampleNames.push_back((string)myParsedHeaderLine[position]);
        position++;
    }
    numSamples=(int)SampleNames.size();
    return true;
}

bool m3vcfHeader::write(IFILE filePtr)
{
    
    if(filePtr == NULL)
    {
        // No file was passed in.
        myStatus.setStatus(StatGenStatus::FAIL_ORDER,
                           "Need to pass in an open file ptr to m3vcfHeader::write.");
        return(false);
    }

    for(std::vector<string>::iterator iter = myHeaderLinesNew.begin();
        iter != myHeaderLinesNew.end(); iter++)
    {
        ifprintf(filePtr, "%s\n", iter->c_str());
        
    }
    return true;
}
// Return the error after a failed call.
const StatGenStatus& m3vcfHeader::getStatus()
{
    return(myStatus);
}


int m3vcfHeader::getNumMetaLines()
{
    int numHeaderLines = myHeaderLinesNew.size();
    if((numHeaderLines >= 1) && (myHasHeaderLine))
    {
        // Remove the header line from the count.
        return(numHeaderLines-1);
    }
    return(numHeaderLines);
}


const char* m3vcfHeader::getMetaLine(unsigned int index)
{
    if(index >= myHeaderLinesNew.size())
    {
        return(NULL);
    }
    else
    {
        return(myHeaderLinesNew[index].c_str());
    }
    return(NULL);
}


const char* m3vcfHeader::getHeaderLine()
{
    // Make sure the last header line is synced with the parsed header line.
    syncHeaderLine();
    if(myHasHeaderLine)
    {
        return(myHeaderLinesNew.back().c_str());
    }
    return(NULL);
}



bool m3vcfHeader::checkMergeHeader(m3vcfHeader &Header)
{
    if(myHeaderLinesNew.size()==0) return true;
    if(numSamples!=Header.numSamples) return false;
    for(int i=0;i<numSamples;i++)
    {
        if(SampleNames[i]!=Header.SampleNames[i]) return false;
//        if(SamplePloidy[i]!=Header.SamplePloidy[i]) return false;
    }
    return true;
}

void m3vcfHeader::mergeHeader(m3vcfHeader &Header)
{
    if(myHeaderLinesNew.size()==0)
    {
        numSamples=Header.numSamples;
        myHeaderLinesNew=Header.myHeaderLinesNew;
        SampleNames=Header.SampleNames;
        SamplePloidy=Header.SamplePloidy;
        myHasHeaderLine=Header.myHasHeaderLine;
    }
    else
    {
        int numHeaderLines = Header.myHeaderLinesNew.size();
        for(int i=0; i<numHeaderLines - 1; i++)
        {
            if(!ifExistsMetaLine(Header.myHeaderLinesNew[i].c_str()))
            {
                appendMetaLine(Header.myHeaderLinesNew[i].c_str());
            }
        }
    }
}

bool m3vcfHeader::ifExistsMetaLine(const char* metaLine)
{
    // Check that the line starts with "##".
    if(strncmp(metaLine, "##", 2) != 0)
    {
        // Does not start with "##"
        return(false);
    }
    if(!myHasHeaderLine)
    {
        // No header line, so just add to the end of the vector.
        return false;
    }

    int numHeaderLines = myHeaderLinesNew.size();
    for(int i=0; i<numHeaderLines; i++)
    {
        if(strcmp(metaLine, myHeaderLinesNew[i].c_str()) == 0)
            return true;
    }
    return false;
}


const char* m3vcfHeader::getSampleName(unsigned int index) const
{
    if(!myHasHeaderLine)
    {
        // No header.
        return(NULL);
    }
    int position = index + NUM_NON_SAMPLE_HEADER_COLS;

    if(position >= myParsedHeaderLine.Length())
    {
        // Out of range.
        return(NULL);
    }

    return(myParsedHeaderLine[position].c_str());
}


int m3vcfHeader::getSampleIndex(const char* sampleName) const
{
    if(!myHasHeaderLine)
    {
        // No header.
        return(-1);
    }
    for(int index = NUM_NON_SAMPLE_HEADER_COLS;
        index < myParsedHeaderLine.Length(); index++)
    {
        if(myParsedHeaderLine[index] == sampleName)
        {
            // Found.
            return(index - NUM_NON_SAMPLE_HEADER_COLS);
        }
    }
    // Not found.
    return(-1);
}



bool m3vcfHeader::appendMetaLine(const char* metaLine)
{
    // Check that the line starts with "##".
    if(strncmp(metaLine, "##", 2) != 0)
    {
        // Does not start with "##"
        return(false);
    }
    if(!myHasHeaderLine)
    {
        // No header line, so just add to the end of the vector.
        myHeaderLinesNew.push_back(metaLine);
        return(true);
    }
    // There is a header line, so insert this just before that line.
    // The headerLine is one position before "end".
    std::vector<string>::iterator headerLine = myHeaderLinesNew.end();
    --headerLine;
    // Insert just before the header line.
    myHeaderLinesNew.insert(headerLine, metaLine);
    return(true);
}


bool m3vcfHeader::addHeaderLine(const char* headerLine)
{
    // Check that the line starts with "#".
    if(strncmp(headerLine, "#", 1) != 0)
    {
        // Does not start with "#"
        return(false);
    }

    if(myHasHeaderLine)
    {
        // There is a header line, so replace the current line.
        myHeaderLinesNew.back() = headerLine;
    }
    else
    {
        // There is not a header line, so add it
        myHeaderLinesNew.push_back(headerLine);
    }

    myHasHeaderLine = true;
    // Parse the header line to get the sample information.
    myParsedHeaderLine.ReplaceColumns(headerLine, '\t');
    return SparseSampleNames();
}


void m3vcfHeader::syncHeaderLine()
{
    if(!myHasHeaderLine)
    {
        // No header line, so nothing to sync.
        return;
    }
    // Get the last header line and see if it is set.
    string& hdrLine = myHeaderLinesNew.back();

    if(hdrLine.length()==0)
    {
        // The header line is not set, so set it.
        for(int i  = 0; i < myParsedHeaderLine.Length(); i++)
        {
            if(i != 0)
            {
                hdrLine += '\t';
            }
            hdrLine += (string)myParsedHeaderLine[i];
        }
    }
}
