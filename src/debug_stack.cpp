// Copyright (c) 2018-2022 California Institute of Technology (“Caltech”) and
// University of Washington. U.S. Government sponsorship acknowledged.
// All rights reserved.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Caltech nor its operating division, the Jet Propulsion
//   Laboratory, nor the names of its contributors may be used to endorse or
//   promote products derived from this software without specific prior written
//   permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifdef _WIN32

void installSegHandler()
{
}

#else

#include <iostream>
#include <execinfo.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ucontext.h>
#include <unistd.h>
#include <cxxabi.h>

using namespace std;


// FROM: https://stackoverflow.com/questions/77005/how-to-automatically-generate-a-stacktrace-when-my-program-crashes

/* This structure mirrors the one found in /usr/include/asm/ucontext.h */
typedef struct _sig_ucontext {
    unsigned long     uc_flags;
    struct ucontext   *uc_link;
    stack_t           uc_stack;
    struct sigcontext uc_mcontext;
    sigset_t          uc_sigmask;
} sig_ucontext_t;

#define TRACE_DEPTH     10

void crit_err_hdlr(int sig, siginfo_t *info, void *ucontext)
{
    sig_ucontext_t *uc = (sig_ucontext_t *)ucontext;

    /* Get the address at the time the signal was raised */
    void *caller_address = 0;
#if defined(__i386__) // gcc specific
    caller_address = (void *)uc->uc_mcontext.eip; // EIP: x86 specific
#elif defined(__x86_64__) // gcc specific
    caller_address = (void *)uc->uc_mcontext.rip; // RIP: x86_64 specific
#else
#error Unsupported architecture. // TODO: Add support for other arch.
#endif

    std::cerr << "signal " << sig 
              << " (" << strsignal(sig) << "), address is " 
              << info->si_addr << " from " << caller_address 
              << std::endl << std::endl;

    void * trace[TRACE_DEPTH];
    int trace_size = backtrace(trace, TRACE_DEPTH);

    trace[1] = caller_address;

    char ** messages = backtrace_symbols(trace, trace_size);    

    printf("[bt] Execution path:\n");
    for(int i = 1; i < trace_size; ++i)
    {
        cerr << "[bt] #" << i << " " << messages[i] << endl;

// This is printing ??:0 for some reason, probably due to addresses being 64bit above (we're using rip not eip cos not i386!)
// So for now, disabling this...
// Also look at: https://stackoverflow.com/questions/3899870/print-call-stack-in-c-or-c/26529030
// demangled = abi::__cxa_demangle(info.dli_sname, NULL, 0, &status);
/*
        // find first occurence of '(' or ' ' in message[i] and assume
        // everything before that is the file name. (Don't go beyond 0 though
        // (string terminator)

        int p = 0;
        while(messages[i][p] != '(' && messages[i][p] != ' ' && messages[i][p] != 0)
        {
            ++p;
        }

        string m = messages[i];
        size_t plus = m.find('+');
        if(plus != string::npos)
        {
            plus++;
        }

        size_t end = m.find(')', plus);
        string addr = m.substr(plus, end-plus);
        cerr << addr << endl;

        //char syscom[256];
        //sprintf(syscom, "addr2line -e %.*s %p", p, messages[i], trace[i]);
        string syscom;
        syscom = "addr2line -e ";
        string loc = messages[i];
        loc = loc.substr(0, p);
        if(loc == "../Piquant")
        {
            loc = "/usr/PIQUANT/Piquant";
        }

        syscom += loc;
        syscom += " ";
        syscom += addr;
        cerr << "calling: " << syscom << endl;
        //last parameter is the file name of the symbol
        int result = system(syscom.c_str());
        if(result != 0)
        {
            cerr << "addr2line result: " << result << endl;
        }
*/
    }

    free(messages);

    exit(EXIT_FAILURE);
}

void installSegHandler()
{
    struct sigaction sigact;

    sigact.sa_sigaction = crit_err_hdlr;
    sigact.sa_flags = SA_RESTART | SA_SIGINFO;

    if (sigaction(SIGSEGV, &sigact, (struct sigaction *)NULL) != 0)
    {
        fprintf(stderr, "error setting signal handler for %d (%s)\n",
        SIGSEGV, strsignal(SIGSEGV));

        exit(EXIT_FAILURE);
    }
}

#endif